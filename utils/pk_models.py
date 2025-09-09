import numpy as np


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
