import numpy as np
from typing import Iterable, List, Tuple, Optional
from .dosing_utils import expand_split_dose


def bateman(t: np.ndarray, dose: float, t0: float, ka: float, ke: float) -> np.ndarray:
    """One-compartment Bateman function for an oral dose starting at t0.

    Returns concentration/amount-like values in the same relative units as `dose`.
    """
    c = np.zeros_like(t, dtype=float)
    m = t >= t0
    tm = t[m] - float(t0)
    # guard for ka ~= ke
    if abs(ka - ke) < 1e-9:
        ka = float(ka) + 1e-6
    c[m] = float(dose) * (ka / (ka - ke)) * (np.exp(-ke * tm) - np.exp(-ka * tm))
    c[c < 0] = 0.0
    return c


def curves_from_schedule(
    t: np.ndarray,
    schedule: Iterable[Tuple[float, float]],
    ka: float,
    ke: float,
) -> List[np.ndarray]:
    """Return Bateman curves for (time, dose) schedule with pre-dose masked.

    Each element in `schedule` is a (t0, dose) pair.
    """
    curves: List[np.ndarray] = [bateman(t, d, td, ka, ke) for td, d in schedule]
    for curve, (td, _) in zip(curves, schedule):
        curve[t < td] = np.nan
    return curves


# === Caffeine model (PK) ===
# Fast absorption vs amphetamine; half-life mid-range ~5 h
KA_CAF: float = 2.0                 # 1/h, adjust to taste
KE_CAF: float = np.log(2) / 5.0     # 1/h, t1/2 ≈ 5 h


def _expand_caffeine_schedule(
    schedule: Optional[Iterable[Tuple[float, float]]]
) -> List[Tuple[float, float, Optional[float], Optional[int]]]:
    """Normalize caffeine schedule entries to a (t0, mg, duration_min?, parts?) form.

    Accepts entries:
    - (t0, mg)
    - (t0, mg, duration_min)
    - (t0, mg, duration_min, parts)
    Returns a list of 4-tuples with duration_min/parts set to None when absent.
    """
    out: List[Tuple[float, float, Optional[float], Optional[int]]] = []
    if not schedule:
        return out
    for item in schedule:
        if not isinstance(item, (tuple, list)):
            continue
        if len(item) == 2:
            t0, mg = float(item[0]), float(item[1])
            out.append((t0, mg, None, None))
        elif len(item) >= 3:
            t0, mg = float(item[0]), float(item[1])
            dur = float(item[2])
            prt = int(item[3]) if len(item) >= 4 and item[3] is not None else None
            out.append((t0, mg, dur, prt))
    return out


def caffeine_total_curve(
    t: np.ndarray,
    schedule: Optional[Iterable[Tuple[float, float]]],
    *,
    ka: float = KA_CAF,
    ke: float = KE_CAF,
) -> Optional[np.ndarray]:
    """Build a cumulative caffeine PK curve from a schedule.

    Each schedule entry can be (t0, mg), (t0, mg, duration_min), or
    (t0, mg, duration_min, parts). Splitting uses 1/min default or `parts`.
    Returns an array matching t, or None if schedule is empty.
    """
    entries = _expand_caffeine_schedule(schedule)
    if not entries:
        return None
    curves: List[np.ndarray] = []
    for t0, mg, dur_min, prt in entries:
        if dur_min is None or dur_min <= 0:
            subs = [(t0, mg)]
        else:
            subs = expand_split_dose(t0, mg, dur_min, parts=prt)
        for ts, d in subs:
            curves.append(bateman(t, d, ts, ka, ke))
    if not curves:
        return None
    return np.nansum(np.vstack([np.nan_to_num(c) for c in curves]), axis=0)
