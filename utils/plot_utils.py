from typing import Tuple


def format_time_12h(h: float) -> str:
    """Format a decimal hour (0–24+) as a 12-hour clock label.

    Examples:
    - 0.0   -> '12am'
    - 12.0  -> '12pm'
    - 13.5  -> '1:30pm'
    - 24.0  -> '12am'
    """
    try:
        hf = float(h)
    except Exception:
        hf = 0.0

    # Round to nearest minute and wrap within 24h for display
    total_min = int(round(hf * 60)) % (24 * 60)
    h24 = total_min // 60
    m = total_min % 60

    suffix = "am" if h24 < 12 else "pm"
    h12 = h24 % 12
    if h12 == 0:
        h12 = 12

    if m == 0:
        return f"{h12}{suffix}"
    return f"{h12}:{m:02d}{suffix}"
