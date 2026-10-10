"""Rules for choosing which patients a two-arm randomised trial should enrol.

All sums run over patient types, weighted by how common each type is among
screened patients.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.stats import norm

_OBJECTIVES = ("power", "cost", "benefit")


@dataclass
class Eligibility:
    """Result of :func:`eligibility`."""

    included: np.ndarray
    ratio: np.ndarray
    threshold: float
    enrolled: float
    mean_effect: float
    ncp: float
    relative_duration: float
    objective: str
    cost_share: Optional[float] = None
    power: Optional[float] = None
    elasticity: Optional[float] = None


def _check(effect, variance, weight, information, prognosis, n_min=1):
    arrays = {"effect": effect, "variance": variance, "weight": weight,
              "information": information, "prognosis": prognosis}
    arrays = {k: None if v is None else np.atleast_1d(np.asarray(v, dtype=float)) for k, v in arrays.items()}
    n = max([n_min] + [a.size for a in arrays.values() if a is not None])
    out = {}
    for name, a in arrays.items():
        if a is None:
            out[name] = None
            continue
        if np.isnan(a).any():
            raise ValueError(f"`{name}` contains missing values.")
        if a.size == 1:
            a = np.repeat(a, n)
        if a.size != n:
            raise ValueError(f"`{name}` must have length 1 or the common length {n}.")
        out[name] = a
    if (out["variance"] <= 0).any():
        raise ValueError("`variance` must be positive.")
    if (out["weight"] < 0).any() or out["weight"].sum() <= 0:
        raise ValueError("`weight` must be non-negative with a positive sum.")
    if out["information"] is not None and (out["information"] <= 0).any():
        raise ValueError("`information` must be positive.")
    out["weight"] = out["weight"] / out["weight"].sum()
    return out


def _signal_noise(effect, variance, information, prognosis, w, included=None, centre=None):
    # Difference in means with 1:1 allocation: signal = effect, noise = 4 * variance,
    # so the non-centrality parameter is sqrt(n_screened) * G / sqrt(V).
    if information is not None:
        return effect * information, information
    v = 4 * variance
    if prognosis is not None and (centre is not None or (included is not None and included.any())):
        mu1 = prognosis + effect
        if centre is None:
            f = w[included].sum()
            centre = ((w * prognosis)[included].sum() / f, (w * mu1)[included].sum() / f)
        v = 4 * (variance + 0.5 * (prognosis - centre[0]) ** 2 + 0.5 * (mu1 - centre[1]) ** 2)
    return effect, v


def _value(G, V, F, H, objective, screen_cost, patient_cost, n_screened, z):
    G, V, F = (np.asarray(a, dtype=float) for a in (G, V, F))
    out = np.full(G.shape, -np.inf)
    ok = (G > 0) & (V > 0)
    if objective == "power":
        out[ok] = np.log(G[ok]) - 0.5 * np.log(V[ok])
    elif objective == "cost":
        out[ok] = -(np.log(V[ok]) - 2 * np.log(G[ok]) + np.log(screen_cost + patient_cost * F[ok]))
    else:
        H = np.asarray(H, dtype=float)
        ok &= H > 0
        nu = np.sqrt(n_screened) * G[ok] / np.sqrt(V[ok])
        out[ok] = norm.logcdf(nu - z) + np.log(H[ok])
    return out


def _best_prefix(index, g, v, w, h, **kw):
    o = np.argsort(-index, kind="stable")
    G, V, F = np.cumsum((w * g)[o]), np.cumsum((w * v)[o]), np.cumsum(w[o])
    H = None if h is None else np.cumsum((w * h)[o])
    val = _value(G, V, F, H, **kw)
    k = int(np.argmax(val))
    inc = np.zeros(g.size, dtype=bool)
    inc[o[: k + 1]] = True
    return inc, val[k]


def _totals(inc, g, v, w, h):
    return ((w * g)[inc].sum(), (w * v)[inc].sum(), w[inc].sum(), None if h is None else (w * h)[inc].sum())


def _set_value(inc, g, v, w, h, **kw):
    if not inc.any():
        return -np.inf
    G, V, F, H = _totals(inc, g, v, w, h)
    return float(_value(G, V, F, H, **kw))


def _marginal(inc, g, v, w, h, objective, screen_cost, patient_cost, n_screened, z):
    G, V, F, H = _totals(inc, g, v, w, h)
    if objective == "power":
        return g / G - v / (2 * V)
    if objective == "cost":
        return 2 * g / G - v / V - patient_cost / (screen_cost + patient_cost * F)
    nu = np.sqrt(n_screened) * G / np.sqrt(V)
    eps = nu * np.exp(norm.logpdf(nu - z) - norm.logcdf(nu - z))
    return eps * (g / G - v / (2 * V)) + h / H


def _solve_fixed(g, v, w, h, objective, screen_cost, patient_cost, n_screened, z):
    kw = dict(objective=objective, screen_cost=screen_cost, patient_cost=patient_cost,
              n_screened=n_screened, z=z)
    pos = w > 0
    if objective == "power" or (objective == "cost" and patient_cost == 0):
        kw["objective"] = "power"
        inc, _ = _best_prefix(np.where(pos, g / v, -np.inf), g, v, w, h, **kw)
        return inc & pos
    # Candidate sets from the two-parameter family implied by the first-order condition.
    if objective == "cost":
        indices = [np.where(pos, (g - a) / v, -np.inf) for a in np.linspace(0, g.max(), 60)]
    else:
        indices = [np.where(pos, (np.cos(t) * g + np.sin(t) * h) / v, -np.inf) for t in np.linspace(0, np.pi / 2, 60)]
    inc, val = None, -np.inf
    for index in indices:
        cand, cval = _best_prefix(index, g, v, w, h, **kw)
        if cval > val:
            inc, val = cand, cval
    inc = inc & pos
    # Polish: apply the rule until the set stops changing, then try single flips.
    for _ in range(200):
        new = (_marginal(inc, g, v, w, h, **kw) >= 0) & pos
        nv = _set_value(new, g, v, w, h, **kw)
        if np.array_equal(new, inc) or not nv > val:
            break
        inc, val = new, nv
    while True:
        G, V, F, H = _totals(inc, g, v, w, h)
        s = np.where(inc, -1.0, 1.0) * w
        vals = _value(G + s * g, V + s * v, F + s, None if h is None else H + s * h, **kw)
        vals[~pos] = -np.inf
        k = int(np.argmax(vals))
        if not vals[k] > val + 1e-12:
            break
        inc[k] = not inc[k]
        val = vals[k]
    return inc


def eligibility(effect, variance=1.0, weight=1.0, information=None, objective="power",
                screen_cost=1.0, patient_cost=0.0, benefit=None, n_screened=None,
                alpha=0.025, prognosis=None) -> Eligibility:
    """Find which patient types a two-arm randomised trial should enrol.

    Parameters
    ----------
    effect : expected treatment effect of each patient type.
    variance : outcome variance within an arm for each type.
    weight : frequency of each type among screened patients.
    information : optional Fisher information per enrolled patient. If given,
        ``effect`` is read as a model coefficient (log odds ratio, log hazard ratio).
    objective : ``"power"`` (highest power for a fixed number screened, or the
        shortest trial for a fixed power), ``"cost"`` (lowest cost for a fixed
        power) or ``"benefit"`` (probability of a positive trial times the total
        benefit of the enrolled types). For ``"power"`` the best set consists of
        the types with the largest effect-to-noise ratios; the function ranks the
        types and compares all candidate sets.
    screen_cost, patient_cost : cost per screened patient (including the cost of
        time) and per enrolled patient.
    benefit : health gain per future patient of each type. Defaults to ``effect``.
    n_screened : number of patients screened during the trial (benefit objective).
    alpha : one-sided significance level.
    prognosis : optional expected outcome under control. If given, the analysis
        is taken to be unadjusted and the set is found by a search that is not
        guaranteed to be exact.

    Examples
    --------
    >>> import numpy as np
    >>> delay = np.arange(0.005, 6, 0.01)
    >>> fit = eligibility(1 - delay / 6)
    >>> round(float(delay[fit.included].max()), 1)
    4.0
    """
    if objective not in _OBJECTIVES:
        raise ValueError(f"`objective` must be one of {_OBJECTIVES}.")
    c = _check(effect, variance, weight, information, prognosis)
    effect, w = c["effect"], c["weight"]
    pos = w > 0
    if not (effect[pos] > 0).any():
        raise ValueError("No patient type has a positive expected effect.")
    if objective == "cost" and (screen_cost < 0 or patient_cost < 0 or screen_cost + patient_cost <= 0):
        raise ValueError("Costs must be non-negative and not both zero.")
    h = None
    if objective == "benefit":
        if n_screened is None or n_screened <= 0:
            raise ValueError("`n_screened` is needed for the benefit objective.")
        h = effect if benefit is None else np.resize(np.asarray(benefit, dtype=float), effect.size)
    z = norm.ppf(1 - alpha)
    profile = c["prognosis"] is not None and c["information"] is None

    def kw(obj):
        return dict(objective=obj, screen_cost=screen_cost, patient_cost=patient_cost,
                    n_screened=n_screened, z=z)

    def solve(obj, included=None, centre=None):
        g, v = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, included, centre)
        return _solve_fixed(g, v, w, h, **kw(obj))

    def value_of(included, obj):
        if not included.any():
            return -np.inf
        g, v = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, included)
        return _set_value(included, g, v, w, h, **kw(obj))

    def find(obj):
        if not profile:
            return solve(obj, included=(effect > 0) & pos)
        # Unadjusted analysis: the noise of a type depends on the mean outcome of
        # the enrolled patients. Search over that mean, then refine.
        start = (effect > 0) & pos
        prog = c["prognosis"][pos]
        if pos.sum() <= 100:
            u = np.unique(prog)
            centres = np.unique(np.concatenate([u, (u[1:] + u[:-1]) / 2]))
        else:
            centres = np.quantile(prog, np.linspace(0.01, 0.99, 65))
        shifts = np.unique(np.quantile(effect[start], [0.1, 0.5, 0.9]))
        starts = [start] + [solve(obj, centre=(m, m + s)) for m in centres for s in shifts]
        best, best_val = None, -np.inf
        seen = set()
        for inc in starts:
            key = inc.tobytes()
            if key in seen:
                continue
            seen.add(key)
            val = value_of(inc, obj)
            while True:
                new = solve(obj, included=inc)
                nv = value_of(new, obj)
                if not nv > val + 1e-12:
                    break
                inc, val = new, nv
            if val > best_val:
                best, best_val = inc.copy(), val
        # Local search: single flips, and swaps when the number of types is small.
        idx = np.flatnonzero(pos)
        improved = True
        while improved:
            improved = False
            for k in idx:
                trial = best.copy()
                trial[k] = not trial[k]
                tv = value_of(trial, obj)
                if tv > best_val + 1e-12:
                    best, best_val, improved = trial, tv, True
            if idx.size <= 60:
                for a in idx[best[idx]]:
                    for b in idx[~best[idx]]:
                        if not best[a] or best[b]:
                            continue
                        trial = best.copy()
                        trial[a], trial[b] = False, True
                        tv = value_of(trial, obj)
                        if tv > best_val + 1e-12:
                            best, best_val, improved = trial, tv, True
        return best

    inc = find(objective)
    ref = inc if objective == "power" else find("power")
    g, v = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, inc)
    G, V, F, _ = _totals(inc, g, v, w, None)
    gr, vr = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, ref)
    Gr, Vr, _, _ = _totals(ref, gr, vr, w, None)
    ratio = g / v
    out = Eligibility(included=inc, ratio=ratio, threshold=float(ratio[inc].min() / (G / V)),
                      enrolled=float(F), mean_effect=float((w * effect)[inc].sum() / F),
                      ncp=float(G / np.sqrt(V)),
                      relative_duration=float((Gr ** 2 / Vr) / (G ** 2 / V)), objective=objective)
    if objective == "cost":
        out.cost_share = float(screen_cost / (screen_cost + patient_cost * F))
    if objective == "benefit":
        nu = np.sqrt(n_screened) * out.ncp
        out.power = float(norm.cdf(nu - z))
        out.elasticity = float(nu * norm.pdf(nu - z) / out.power)
    return out


def _included(included, n, w):
    if included is None:
        included = np.ones(n, dtype=bool)
    included = np.atleast_1d(np.asarray(included, dtype=bool))
    if included.size != n:
        raise ValueError("`included` must have one value per patient type.")
    return included & (w > 0)


def design_power(effect, variance=1.0, weight=1.0, information=None, prognosis=None,
                 included=None, *, n_screened, alpha=0.025) -> dict:
    """Power of a trial with a given eligibility set (one-sided test)."""
    n_min = 1 if included is None else np.size(included)
    c = _check(effect, variance, weight, information, prognosis, n_min)
    effect, w = c["effect"], c["weight"]
    inc = _included(included, effect.size, w)
    if not inc.any():
        raise ValueError("No patient type is enrolled.")
    g, v = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, inc)
    G, V, F, _ = _totals(inc, g, v, w, None)
    ncp = np.sqrt(n_screened) * G / np.sqrt(V)
    return {"enrolled": float(F), "mean_effect": float((w * effect)[inc].sum() / F),
            "ncp": float(ncp), "power": float(norm.cdf(ncp - norm.ppf(1 - alpha)))}


def relative_duration(effect, variance=1.0, weight=1.0, information=None, prognosis=None,
                      included=None) -> float:
    """Duration of a trial with a given eligibility set, relative to the power-optimal set."""
    n_min = 1 if included is None else np.size(included)
    c = _check(effect, variance, weight, information, prognosis, n_min)
    effect, w = c["effect"], c["weight"]
    inc = _included(included, effect.size, w)
    best = eligibility(effect, c["variance"], w, c["information"], prognosis=c["prognosis"])
    g, v = _signal_noise(effect, c["variance"], c["information"], c["prognosis"], w, inc)
    G, V, _, _ = _totals(inc, g, v, w, None)
    if G <= 0:
        return float("inf")
    return float(best.ncp ** 2 / (G ** 2 / V))


def inclusion_bound(effect, variance=1.0, weight=1.0, information=None) -> float:
    """Upper bound ``1 + CV**2`` on the relative duration of a trial that enrols everyone."""
    c = _check(effect, variance, weight, information, None)
    g, v = _signal_noise(c["effect"], c["variance"], c["information"], None, c["weight"])
    m = c["weight"] * v
    r = g / v
    m1 = (m * r).sum() / m.sum()
    if m1 <= 0:
        return float("inf")
    return float((m * r ** 2).sum() / m.sum() / m1 ** 2)


def power_elasticity(power=0.8, alpha=0.025):
    """Percentage change in power per one percent change in the non-centrality parameter."""
    power = np.asarray(power, dtype=float)
    if ((power <= 0) | (power >= 1)).any():
        raise ValueError("`power` must lie strictly between 0 and 1.")
    zb = norm.ppf(power)
    out = (norm.ppf(1 - alpha) + zb) * norm.pdf(zb) / power
    return float(out) if out.ndim == 0 else out


def rule_threshold(objective="power", power=0.8, alpha=0.025, cost_share=1.0):
    """Effect of the marginal type at the best set, as a fraction of the average effect among enrolled patients.

    ``"power"`` and ``"cost"`` give ``1 - cost_share / 2``; ``"benefit"`` gives
    ``e / (2 * (1 + e))`` and ``"coverage"`` gives ``1/2 - 1/e``, where ``e`` is
    the elasticity of power. All types are assumed to be equally noisy and small;
    for a whole group use :func:`subgroup_threshold`.
    """
    if objective in ("power", "cost"):
        cs = np.asarray(cost_share, dtype=float)
        if ((cs < 0) | (cs > 1)).any():
            raise ValueError("`cost_share` must lie between 0 and 1.")
        out = 1 - cs / 2
        return float(out) if out.ndim == 0 else out
    if objective not in ("benefit", "coverage"):
        raise ValueError("`objective` must be 'power', 'cost', 'benefit' or 'coverage'.")
    e = power_elasticity(power, alpha)
    return e / (2 * (1 + e)) if objective == "benefit" else 0.5 - 1 / e


def subgroup_threshold(rho, reference="others"):
    """Exact threshold for adding or dropping a whole group of patients.

    ``rho`` is the information of the group divided by that of the other enrolled
    patients. With ``reference="others"`` the result is ``1 / (1 + sqrt(1 + rho))``:
    adding the group shortens the trial if and only if its effect-to-noise ratio,
    divided by the pooled ratio of the other enrolled patients, is above it. With
    ``reference="all"`` the result is ``sqrt(1 + rho) / (1 + sqrt(1 + rho))``: an
    enrolled group should be dropped if and only if its ratio, divided by the pooled
    ratio of all enrolled patients including the group, is below it.
    """
    if reference not in ("others", "all"):
        raise ValueError("`reference` must be 'others' or 'all'.")
    rho = np.asarray(rho, dtype=float)
    if (rho < 0).any():
        raise ValueError("`rho` must be non-negative.")
    s = np.sqrt(1 + rho)
    out = 1 / (1 + s) if reference == "others" else s / (1 + s)
    return float(out) if out.ndim == 0 else out


def time_window(effect: Callable, density: Optional[Callable] = None, *, upper: float,
                lower: float = 0.0, variance=1.0, n_grid: int = 4000, **kwargs) -> dict:
    """Best treatment time window when the effect depends on the delay to treatment.

    ``effect``, ``density`` and (optionally) ``variance`` are functions of the delay.
    Further keyword arguments are passed to :func:`eligibility`.
    """
    step = (upper - lower) / n_grid
    t = lower + (np.arange(n_grid) + 0.5) * step
    v = variance(t) if callable(variance) else variance
    wt = np.ones_like(t) if density is None else density(t)
    fit = eligibility(effect(t), variance=v, weight=wt, **kwargs)
    return {"delay": t, "fit": fit, "close": float(t[fit.included].max() + step / 2)}
