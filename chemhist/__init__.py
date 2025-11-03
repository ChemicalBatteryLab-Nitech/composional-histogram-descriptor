# chemhist/__init__.py

from importlib import resources
import os

from .core import get_descriptor

__all__ = ["get_descriptor"]

# --- Safe data directory detection for both normal and Store Python ---
try:
    # Preferred modern API
    DATA_PATH = resources.files("chemhist") / "data"
except Exception:
    # Fallback: use relative directory or current working dir
    here = os.path.abspath(os.path.dirname(__file__) if "__file__" in globals() else os.getcwd())
    DATA_PATH = os.path.join(here, "data")

# Export the path (for developers)
__data_path__ = str(DATA_PATH)
