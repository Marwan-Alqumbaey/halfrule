import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.stats import norm

from halfrule import (eligibility, inclusion_bound, power_elasticity, relative_duration,
                      rule_threshold, subgroup_threshold, time_window)


def test_thresholds_match_closed_forms():
    assert rule_threshold("power") == 0.5
    assert rule_threshold("cost", cost_share=0.5) == 0.75
    e = power_elasticity(0.8)
    assert e == pytest.approx((norm.ppf(0.975) + norm.ppf(0.8)) * norm.pdf(norm.ppf(0.8)) / 0.8)
    assert rule_threshold("benefit", power=0.8) == pytest.approx(e / (2 * (1 + e)))
    assert round(rule_threshold("benefit", power=0.8), 3) == 0.248
    assert rule_threshold("coverage", power=0.8) < 0
    with pytest.raises(ValueError):
        rule_threshold("cost", cost_share=2)
    with pytest.raises(ValueError):
        power_elasticity(1)


def test_subgroup_threshold():
    assert subgroup_threshold(0) == 0.5
    prev = np.array([0.25, 0.5, 0.75])
    assert np.allclose(subgroup_threshold((1 - prev) / prev), np.sqrt(prev) / (1 + np.sqrt(prev)))
    with pytest.raises(ValueError):
        subgroup_threshold(-1)
    with pytest.raises(ValueError):
        subgroup_threshold(1, reference="best")


def test_group_thresholds_agree_with_direct_comparison():
    rng = np.random.default_rng(1)
    kappa = rng.uniform(-0.5, 1.5, 2000)
    rho = np.exp(rng.uniform(np.log(0.01), np.log(20), 2000))
    better = (1 + kappa * rho) / np.sqrt(1 + rho) > 1
    assert np.array_equal(better, kappa > subgroup_threshold(rho))
    keep = 1 + kappa * rho > 0
    kappa_all = kappa * (1 + rho) / (1 + kappa * rho)
    assert np.array_equal(better[keep], kappa_all[keep] > subgroup_threshold(rho[keep], reference="all"))
    assert (subgroup_threshold(rho) <= 0.5).all()
    assert (subgroup_threshold(rho, reference="all") >= 0.5).all()


def test_linear_decay_closes_at_two_thirds():
    assert time_window(lambda t: 1 - t / 6.3, upper=8)["close"] == pytest.approx(4.2, abs=2e-3)
    assert time_window(lambda t: 1 - t / 9, upper=9)["close"] == pytest.approx(6, abs=2e-3)


def test_exponential_decay():
    u = brentq(lambda x: np.exp(x) - 1 - 2 * x, 0.5, 3)
    assert time_window(lambda t: np.exp(-t / 2), upper=20)["close"] == pytest.approx(2 * u, abs=5e-3)


def test_benefit_objective_gives_quarter_rule():
    t = np.arange(0.0005, 1, 0.001)
    fit = eligibility(1 - t, objective="benefit", n_screened=60)
    e = fit.elasticity
    assert fit.threshold == pytest.approx(e / (2 * (1 + e)), abs=5e-3)
    assert fit.enrolled > eligibility(1 - t).enrolled


def test_cost_objective_gives_cost_rule():
    t = np.arange(0.0005, 1, 0.001)
    fit = eligibility(1 - t, objective="cost", screen_cost=1, patient_cost=3)
    assert fit.threshold == pytest.approx(1 - fit.cost_share / 2, abs=5e-3)
    free = eligibility(1 - t, objective="cost", screen_cost=1, patient_cost=0)
    assert np.array_equal(free.included, eligibility(1 - t).included)


def test_constant_effect_follows_information():
    rng = np.random.default_rng(1)
    info = rng.gamma(2, 0.25, 4000)
    fit = eligibility(0.3, information=info, objective="cost", screen_cost=1, patient_cost=5)
    cut = (1 - fit.cost_share) * info[fit.included].mean()
    assert info[fit.included].min() == pytest.approx(cut, rel=0.02)
    assert eligibility(0.3, information=info).included.all()


def test_two_subgroup_table():
    # Simon and Maitournam (2004): screened patients, untargeted relative to targeted design.
    for prev, published in ((0.75, 0.98), (0.5, 0.89), (0.25, 0.64)):
        ratio = (relative_duration([1, 0.5], weight=[prev, 1 - prev])
                 / relative_duration([1, 0.5], weight=[prev, 1 - prev], included=[True, False]))
        assert round(ratio, 2) == published
    assert eligibility([1, 0.34], weight=[0.25, 0.75]).included.all()
    assert eligibility([1, 0.32], weight=[0.25, 0.75]).included.tolist() == [True, False]


def test_inclusion_bound():
    rng = np.random.default_rng(14)
    for _ in range(20):
        g, v, w = rng.exponential(size=200), rng.uniform(0.5, 2, 200), rng.uniform(size=200)
        assert relative_duration(g, v, w) <= inclusion_bound(g, v, w) + 1e-12
    t = np.arange(0.0005, 1, 0.001)
    assert inclusion_bound(1 - t) == pytest.approx(4 / 3, abs=1e-4)
