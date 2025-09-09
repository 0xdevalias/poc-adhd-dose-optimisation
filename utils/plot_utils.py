from typing import Tuple


def label_hour(h: float) -> str:
    """Format a decimal hour (0–24+) as a 12-hour clock label like '8am'."""
    try:
        h24 = int(float(h) % 24)
    except Exception:
        h24 = 0
    suffix = "am" if h24 < 12 else "pm"
    h12 = h24 % 12
    if h12 == 0:
        h12 = 12
    return f"{h12}{suffix}"

