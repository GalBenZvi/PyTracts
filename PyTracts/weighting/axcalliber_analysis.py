from PyTracts.weighting.relevant_paths import RelevantPaths
import matlab.engine
from pathlib import Path

MATLAB_CODES_DIR = RelevantPaths.matlab_dir.value


class AxCalliberAnalysis:
    def __init__(
        self,
        dwi_fname: Path,
        bval_fname: Path,
        bvec_fname: Path,
        matlab_dir: Path = MATLAB_CODES_DIR,
        small_delta: float = 15.5,
        big_delta: float = 60,
        g_max: float = 7.2,
    ):
        self.dwi_fname = dwi_fname
        self.bval_fname = bval_fname
        self.bvec_fname = bvec_fname
        self.matlab_dir = matlab_dir
        self.small_delta = float(small_delta)
        self.big_delta = float(big_delta)
        self.g_max = float(g_max)

    def init_matlab_engine(self):
        print("Initiating MATLAB engine...")
        eng = matlab.engine.start_matlab()
        print(self.matlab_dir)
        eng.addpath(str(self.matlab_dir),nargout=0)
        return eng

    def perform_analysis(self, eng):
        check = list(Path(self.dwi_fname.parent).glob("*1.5_2_AxPasi5*"))
        if not check:
            print("Calling AxCaliberV6 function with inputs:")
            print(f"dwi_fname = {self.dwi_fname}")
            print(f"bval_fname = {self.bval_fname}")
            print(f"bvec_fname = {self.bvec_fname}")
            print("This may take a while, so grab a cup of coffee...")
            eng.AxCaliberV6(
                str(self.dwi_fname),
                str(self.bval_fname),
                str(self.bvec_fname),
                self.small_delta,
                self.big_delta,
                self.g_max,
                nargout=0,
            )

    def run(self):
        eng = self.init_matlab_engine()
        self.perform_analysis(eng)
