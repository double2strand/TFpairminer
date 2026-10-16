"""Alias of mining/catlas_pipeline.py (shared motif, matching and FDR helpers). No CATlas data or results are used."""
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent.parent))
from catlas_pipeline import *  # noqa: F401,F403
from catlas_pipeline import command_compile  # noqa: F401
