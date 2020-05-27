from PyTracts.weighting.relevant_paths import RelevantPaths
from PyTracts.weighting.axcalliber_analysis import AxCalliberAnalysis
from pathlib import Path
from PyTracts.utils import FSLOUTTYPE
import glob


class AxCalliberTractography:
    def __init__(
        self,
        derivatives_dir: Path,
        subj: str = None,
        weight_by: str = "1.5_2_AxPasi5",
        preprocessed_fname: str = "acq-AP_dwi_preprocessed_biascorr",
        small_delta: float = 15.5,
        big_delta: float = 60,
        g_max: float = 7.2,
    ):
        self.derivatives = Path(derivatives_dir)
        if subj:
            subjects = [subj]
        else:
            subjects = [subj.name for subj in self.derivatives.glob("sub-*")]
        subjects.sort()
        subjects_dict = dict()
        for subj in subjects:
            subjects_dict[subj] = self.derivatives / subj
        self.subjects = subjects_dict
        self.weight_vy = weight_by
        self.preprocessed_fname = preprocessed_fname
        self.small_delta = small_delta
        self.big_delta = big_delta
        self.g_max = g_max

    def init_subject_params(self, folder_name: Path):
        dwi_folder = folder_name / "dwi"
        dwi_fname = list(dwi_folder.glob(f"*{self.preprocessed_fname}{FSLOUTTYPE}"))[0]
        bvec_fname = list(dwi_folder.glob("*.bvec"))[0]
        bval_fname = list(dwi_folder.glob("*.bval"))[0]
        return dwi_fname, bval_fname, bvec_fname

    def perform_axcalliber(self, dwi_fname: Path, bval_fname: Path, bvec_fname: Path):
        axcalliber_analysis = AxCalliberAnalysis(
            dwi_fname,
            bval_fname,
            bvec_fname,
            small_delta=self.small_delta,
            big_delta=self.big_delta,
            g_max=self.g_max,
        )
        axcalliber_analysis.run()

    def run(self):
        for subj in self.subjects:
            print(f"Working on {subj}...")
            folder_name = self.subjects[subj]
            dwi_fname, bval_fname, bvec_fname = self.init_subject_params(folder_name)
            self.perform_axcalliber(dwi_fname, bval_fname, bvec_fname)
