from pathlib import Path
from PyTracts.DipyTractography import (
    dipy_methods,
    GenerateGrayAndWhite,
    LoadData,
    CreateModel,
    CalculateSeeds,
)
from dipy.io.image import load_nifti, load_nifti_data
from dipy.core.gradients import GradientTable
import numpy as np


class GenerateTractsDipy:
    def __init__(
        self,
        mother_dir: Path,
        reconstruction: str = "probabilistic",
        subj: str = None,
        white_label: int = 3,
        gray_label: int = 2,
        sh_order: int = 6,
        roi_radius: int = 10,
        fa_thr: float = 0.7,
        relative_peak_threshold: float = 0.8,
        min_separation_angle: int = 45,
        stopping_threshold: float = 0.18,
        seeds_density: list = [2, 2, 2],
    ):
        self.seeds_density = seeds_density
        self.stopping_threshold = stopping_threshold
        self.sh_order = sh_order
        self.roi_radius = roi_radius
        self.fa_thr = fa_thr
        self.relative_peak_threshold = relative_peak_threshold
        self.min_separation_angle = min_separation_angle
        self.reconstruction = reconstruction
        self.mother_dir = mother_dir
        if subj:
            subjects = [subj]
        else:
            subjects = [subj.name for subj in mother_dir.glob("sub-*")]
        subjects.sort()
        subjects_dict = dict()
        for subj in subjects:
            subjects_dict[subj] = mother_dir / subj
        self.subjects = subjects_dict
        self.white_label = white_label
        self.gray_label = gray_label

    def load_file_names(self, folder_name: Path):
        dwi_folder = folder_name / "dwi"
        for file in dwi_folder.iterdir():
            if file.suffix == ".bvec":
                bvec_file = file
            elif "preprocessed_biascorr.nii.gz" in str(file):
                dwi_file = file
        segment_file_name = folder_name / "anat" / "prep.anat" / "T1_fast_seg.nii.gz"
        bval_file = bvec_file.with_suffix(".bval")
        return dwi_file, segment_file_name, bvec_file, bval_file

    def white_and_gray_masks(self, segment_file_name: Path, dwi_file: Path):
        masks_generator = GenerateGrayAndWhite(
            segment_file_name, dwi_file, self.white_label, self.gray_label
        )
        white_mask, gray_mask = masks_generator.run()
        return white_mask, gray_mask

    def load_data(self, dwi_file: Path, bvec_file: Path, bval_file: Path):
        data_loader = LoadData(dwi_file, bvec_file, bval_file)
        data, affine, hardi_img, gtab = data_loader.run()
        return data, affine, hardi_img, gtab

    def create_model(
        self,
        folder_name: Path,
        data: np.ndarray,
        white_mask: np.ndarray,
        gtab: GradientTable,
    ):
        model_generator = CreateModel(
            folder_name,
            self.reconstruction,
            data,
            gtab,
            white_mask,
            self.sh_order,
            self.roi_radius,
            self.fa_thr,
            self.relative_peak_threshold,
            self.min_separation_angle,
            self.stopping_threshold,
        )
        model = model_generator.run()
        return model

    def calculate_seeds(self, gray_mask: np.ndarray, affine: np.ndarray):
        seeds_calculator = CalculateSeeds(gray_mask, affine, self.seeds_density)
        seeds = seeds_calculator.calculate_seeds()
        return seeds

