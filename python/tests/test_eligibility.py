import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm

from halfrule import design_power, eligibility, inclusion_bound, relative_duration

FIXTURES = json.loads((Path(__file__).parent / "fixtures.json").read_text())


def brute(value, k):
    best, arg = -np.inf, None
    for bits in itertools.product([False, True], repeat=k):
        inc = np.array(bits)
        if not inc.any():
            continue
        val = value(inc)
        if val > best:
            best, arg = val, inc
    return best, arg


def test_power_equals_exhaustive_search():
    rng = np.random.default_rng(11)
    for _ in range(15):
        k = 10
        g, v, w = rng.normal(0.3, 0.6, k), rng.uniform(0.3, 3, k), rng.uniform(size=k)
        if not (g > 0).any():
            continue
        fit = eligibility(g, v, w)
        val = lambda inc: (w * g)[inc].sum() / np.sqrt((w * v)[inc].sum())
        assert val(fit.included) == pytest.approx(brute(val, k)[0])
        assert fit.threshold >= 0.5


def test_cost_and_benefit_equal_exhaustive_search():
    rng = np.random.default_rng(12)
    for _ in range(15):
        k = 10
        g, v, w, h = rng.normal(0.3, 0.6, k), rng.uniform(0.3, 3, k), rng.uniform(size=k), rng.uniform(0, 2, k)
        w = w / w.sum()
        if not (g > 0).any():
            continue

        def cost(inc):
            G = (w * g)[inc].sum()
            return -np.inf if G <= 0 else -(4 * (w * v)[inc].sum() / G ** 2) * (1 + 2 * w[inc].sum())

        fit = eligibility(g, v, w, objective="cost", screen_cost=1, patient_cost=2)
        assert cost(fit.included) == pytest.approx(brute(cost, k)[0])

        def ben(inc):
            G = (w * g)[inc].sum()
            if G <= 0:
                return -np.inf
            return norm.cdf(np.sqrt(50) * G / np.sqrt(4 * (w * v)[inc].sum()) - norm.ppf(0.975)) * (w * h)[inc].sum()

        fit = eligibility(g, v, w, objective="benefit", benefit=h, n_screened=50)
        assert ben(fit.included) == pytest.approx(brute(ben, k)[0])


def test_unadjusted_analysis():
    rng = np.random.default_rng(13)
    for _ in range(15):
        k = 10
        g, v, w, mu = np.abs(rng.normal(0.3, 0.3, k)), rng.uniform(0.5, 2, k), rng.uniform(size=k), rng.normal(size=k)

        def val(inc):
            f = w[inc].sum()
            m0, m1 = (w * mu)[inc].sum() / f, (w * (mu + g))[inc].sum() / f
            noise = v + 0.5 * (mu - m0) ** 2 + 0.5 * (mu + g - m1) ** 2
            return (w * g)[inc].sum() / np.sqrt((w * noise)[inc].sum())

        fit = eligibility(g, v, w, prognosis=mu)
        assert val(fit.included) == pytest.approx(brute(val, k)[0], rel=1e-8)
    x = np.linspace(-3, 3, 601)
    assert eligibility(1, 1, norm.pdf(x)).included.all()
    assert eligibility(1, 1, norm.pdf(x), prognosis=2 * x).enrolled < 1


@pytest.mark.parametrize("case", FIXTURES)
def test_agrees_with_r_package(case):
    g, v, w = (np.array(case[k]) for k in ("effect", "variance", "weight"))
    a = eligibility(g, v, w)
    assert np.flatnonzero(a.included).tolist() == case["power"]["included"]
    assert a.ncp == pytest.approx(case["power"]["ncp"])
    assert a.threshold == pytest.approx(case["power"]["threshold"])
    b = eligibility(g, v, w, objective="cost", screen_cost=1, patient_cost=2.5)
    assert np.flatnonzero(b.included).tolist() == case["cost"]["included"]
    assert b.cost_share == pytest.approx(case["cost"]["cost_share"])
    assert b.relative_duration == pytest.approx(case["cost"]["relative_duration"])
    c = eligibility(g, v, w, objective="benefit", benefit=case["benefit"], n_screened=80)
    assert np.flatnonzero(c.included).tolist() == case["benefit_fit"]["included"]
    assert c.power == pytest.approx(case["benefit_fit"]["power"])
    assert c.elasticity == pytest.approx(case["benefit_fit"]["elasticity"])
    d = eligibility(np.abs(g), v, w, prognosis=case["prognosis"])
    assert np.flatnonzero(d.included).tolist() == case["unadjusted"]["included"]
    assert d.ncp == pytest.approx(case["unadjusted"]["ncp"])
    e = eligibility(g, information=case["information"], weight=w, objective="cost", screen_cost=1, patient_cost=4)
    assert np.flatnonzero(e.included).tolist() == case["score_cost"]["included"]
    assert relative_duration(g, v, w) == pytest.approx(case["all_duration"])
    assert inclusion_bound(g, v, w) == pytest.approx(case["bound"])
    assert design_power(g, v, w, included=g > 0.2, n_screened=300)["power"] == pytest.approx(case["design"])


def test_design_power_matches_textbook_formula():
    out = design_power(0.3, 1, n_screened=400)
    assert out["ncp"] == pytest.approx(0.3 * np.sqrt(400) / 2)
    assert out["power"] == pytest.approx(norm.cdf(3 - norm.ppf(0.975)))
    with pytest.raises(ValueError):
        design_power([1, 2], included=[False, False], n_screened=10)


def test_inputs_are_validated():
    for bad in (dict(effect=[-1, -2]), dict(effect=[1, np.nan]), dict(effect=[1, 2], variance=[1, -1]),
                dict(effect=[1, 2], weight=[1, 2, 3]), dict(effect=[1, 2], objective="benefit"),
                dict(effect=[1, 2], objective="other")):
        with pytest.raises(ValueError):
            eligibility(**bad)
    with pytest.raises(ValueError):
        relative_duration([1, 2], included=[True, False, True])
    assert relative_duration(1, weight=[1, 1], included=[True, False]) == pytest.approx(2)
