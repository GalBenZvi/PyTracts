from pathlib import Path
from enum import Enum

CURRENT_DIR = Path(__file__).parent


class RelevantPaths(Enum):
    matlab_dir = CURRENT_DIR / "matlab_codes"
