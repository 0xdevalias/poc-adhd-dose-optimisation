import sys
from pathlib import Path
from typing import Optional, Dict, Any


def _prompt_yes_no(question: str, default: bool = False) -> bool:
    """Prompt user for a yes/no answer; default used if not a TTY.

    Returns True only for explicit 'y'/'yes' (case-insensitive).
    """
    if not sys.stdin.isatty():
        return default
    try:
        resp = input(f"{question} [y/N]: ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        return default
    return resp in {"y", "yes"}


def ensure_save_path(path: Optional[str], *, on_exist: str = 'prompt', mkdirs: bool = False) -> Optional[Path]:
    """Safety checks for saving to a file path.

    - on_exist: 'prompt' (ask before overwrite), 'overwrite' (always), 'error' (raise if exists)
    - mkdirs: create parent directories if needed

    Returns resolved Path if it is safe to write, or None if user declined (prompt case).
    Raises FileNotFoundError/IsADirectoryError/FileExistsError as appropriate.
    """
    if not path:
        return None

    p = Path(path).expanduser()
    parent = p.parent

    if not parent.exists():
        if mkdirs:
            parent.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Directory does not exist: {parent}")

    if not parent.is_dir():
        raise NotADirectoryError(f"Parent is not a directory: {parent}")

    if p.exists() and p.is_dir():
        raise IsADirectoryError(f"Refusing to overwrite directory: {p}")

    if p.exists() and p.is_file():
        if on_exist == 'overwrite':
            return p

        if on_exist == 'error':
            raise FileExistsError(f"File already exists: {p}")

        # 'prompt'
        # In non-interactive contexts, be strict and raise for safety
        if not sys.stdin.isatty():
            raise FileExistsError(f"File already exists (non-interactive): {p}")

        if _prompt_yes_no(f"File exists: {p}. Overwrite?", default=False):
            return p

        # declined
        return None

    return p


def save_figure_safely(fig, path: Optional[str], *, on_exist: str = 'prompt', mkdirs: bool = False, default_kwargs: Optional[Dict[str, Any]] = None, **savefig_kwargs) -> bool:
    """Save a Matplotlib figure with safety checks and savefig-style kwargs.

    - Safety: validates directory, handles existing files based on `on_exist`.
    - Defaults: applies preferred savefig defaults (e.g., bbox_inches='tight'), caller can override.
    - Flexibility: forwards arbitrary savefig kwargs (format, dpi, metadata, etc.).

    Returns True if saved, False if skipped (e.g., user declined overwrite or no path).
    In non-interactive mode and on_exist='prompt', raises FileExistsError to be safe.
    """
    p = ensure_save_path(path, on_exist=on_exist, mkdirs=mkdirs)
    if not p:
        return False

    # Preferred defaults
    final_kwargs: Dict[str, Any] = {"bbox_inches": "tight"}

    if default_kwargs:
        final_kwargs.update(default_kwargs)

    # Allow caller overrides
    final_kwargs.update(savefig_kwargs)

    fig.savefig(str(p), **final_kwargs)
    return True
