from pathlib import Path
from PyTracts.DipyTractography import (
    dipy_methods,
    GenerateGrayAndWhite,
    LoadData,
    CreateModel,
    CalculateSeeds,
    StreamlineGenerator,
    DefaultStreamlines,
)
from dipy.io.image import load_nifti, load_nifti_data
from dipy.core.gradients import GradientTable
import numpy as np
import glob


class GenerateTractsDipy:
    def __init__(
        self,
        mother_dir: Path,
        bids_dir: Path = None,
        reconstruction: str = "deterministic",
        subj: str = None,
        white_label: int = 3,
        gray_label: int = 2,
        sh_order: int = 6,
        roi_radius: int = 10,
        fa_thr: float = 0.7,
        relative_peak_threshold: float = 0.8,
        min_separation_angle: int = 45,
        stopping_threshold: float = 0.15,
        seeds_density=1,
        tractogram_fname: str = "Whole_brain_tractography.trk",
        sphere: str = "default",
        max_angle: float = 30.0,
        step_size: float = 1,
        small_delta: float = 15.5,
    ):
        self.tractogram_fname = tractogram_fname
        self.sphere = sphere
        self.small_delta = small_delta
        self.max_angle = max_angle
        self.step_size = step_size
        self.seeds_density = seeds_density
        self.stopping_threshold = stopping_threshold
        self.sh_order = sh_order
        self.roi_radius = roi_radius
        self.fa_thr = fa_thr
        self.relative_peak_threshold = relative_peak_threshold
        self.min_separation_angle = min_separation_angle
        self.reconstruction = reconstruction
        self.mother_dir = Path(mother_dir)
        if not bids_dir:
            self.bids_dir = self.mother_dir.parent / "bids_dataset"
        if subj:
            subjects = [subj]
        else:
            subjects = [subj.name for subj in self.mother_dir.glob("sub-*")]
        subjects.sort()
        subjects_dict = dict()
        for subj in subjects:
            subjects_dict[subj] = self.mother_dir / subj
        self.subjects = subjects_dict
        self.white_label = white_label
        self.gray_label = gray_label

    def load_file_names(self, folder_name: Path, subject_bids: Path):
        dwi_folder = folder_name / "dwi"
        fmap_folder = folder_name / "fmap"
        for file in dwi_folder.iterdir():
            if file.suffix == ".bvec":
                bvec_file = file
            elif "preprocessed_biascorr.nii.gz" in str(file):
                dwi_file = file
        for file in fmap_folder.iterdir():
            if "fieldmap_magnitude.nii.gz" in str(file):
                fieldmap_mag = file
            elif "fieldmap_rad" in str(file):
                fieldmap = file
        phasediff_json = glob.glob(f"{str(subject_bids / 'fmap')}/*.json")[0]
        wm_fname = folder_name / "anat" / "prep.anat" / "T1_fast_pve_2.nii.gz"
        gm_fname = folder_name / "anat" / "prep.anat" / "T1_fast_pve_1.nii.gz"
        structual_fname = folder_name / "anat" / "prep.anat" / "T1_biascorr.nii.gz"
        bval_file = bvec_file.with_suffix(".bval")
        return (
            structual_fname,
            dwi_file,
            wm_fname,
            gm_fname,
            bvec_file,
            bval_file,
            fieldmap,
            fieldmap_mag,
            phasediff_json,
        )

    def white_and_gray_masks(
        self,
        structural_fname: Path,
        wm_fname: Path,
        gm_fname: Path,
        dwi_fname: Path,
        fieldmap: Path,
        fieldmap_magnitude: Path,
        phasediff_json: Path,
    ):
        masks_generator = GenerateGrayAndWhite(
            structural_fname,
            wm_fname,
            gm_fname,
            dwi_fname,
            fieldmap,
            fieldmap_magnitude,
            phasediff_json,
        )
        white_mask, gray_mask = masks_generator.run()
        return white_mask, gray_mask

    def load_data(self, dwi_file: Path, bvec_file: Path, bval_file: Path):
        data_loader = LoadData(dwi_file, bvec_file, bval_file, self.small_delta)
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
        csd_fit, stopping_criterion, csa_model = model_generator.run()
        return csd_fit, stopping_criterion, csa_model

    def calculate_seeds(self, gray_mask: np.ndarray, affine: np.ndarray):
        seeds_calculator = CalculateSeeds(gray_mask, affine, self.seeds_density)
        seeds = seeds_calculator.calculate_seeds()
        return seeds

    def generate_streamlines(
        self,
        folder_name: Path,
        csd_fit,
        hardi_img,
        stopping_criterion,
        seeds,
        affine,
        csa_peaks,
    ):
        streamline_generator = StreamlineGenerator(
            folder_name,
            csd_fit,
            hardi_img,
            self.reconstruction,
            stopping_criterion,
            seeds,
            affine,
            self.tractogram_fname,
            self.sphere,
            self.max_angle,
            self.step_size,
            csa_peaks,
        )
        streamlines = streamline_generator.run()
        return streamlines

    def default_streamlines(
        self, folder_name: Path, csa_model, data: np.ndarray, white_matter: np.ndarray
    ):
        csa_calculator = DefaultStreamlines(
            folder_name,
            csa_model,
            data,
            white_matter,
            self.relative_peak_threshold,
            self.min_separation_angle,
        )
        csa_peaks = csa_calculator.generate_peaks()
        return csa_peaks

    def run(self):
        for subj in self.subjects:
            folder_name = self.subjects[subj]
            subject_bids = self.bids_dir / subj
            tracts_dir = folder_name / "tractography"
            if not tracts_dir.is_dir():
                tracts_dir.mkdir()
            print("loading files")
            (
                structual_fname,
                dwi_file,
                wm_fname,
                gm_fname,
                bvec_file,
                bval_file,
                fieldmap,
                fieldmap_mag,
                phasediff_json,
            ) = self.load_file_names(folder_name, subject_bids)
            print("epi_reg + segmentation")
            white_mask, gray_mask = self.white_and_gray_masks(
                structual_fname,
                wm_fname,
                gm_fname,
                dwi_file,
                fieldmap,
                fieldmap_mag,
                phasediff_json,
            )
            print("loading data")
            data, affine, hardi_img, gtab = self.load_data(
                dwi_file, bvec_file, bval_file
            )
            print("generating model")
            csd_fit, stopping_criterion, csa_model = self.create_model(
                folder_name / "tractography", data, white_mask, gtab
            )
            print("calculating seeds")
            seeds = self.calculate_seeds(gray_mask, affine)
            if self.reconstruction.lower() == "default":
                csa_peaks = self.default_streamlines(
                    folder_name / "tractography", csa_model, data, white_mask
                )
            else:
                csa_peaks = None
            print("generating streamlines")
            streamlines = self.generate_streamlines(
                folder_name / "tractography",
                csd_fit,
                hardi_img,
                stopping_criterion,
                seeds,
                affine,
                csa_peaks,
            )


if __name__ == "__main__":
    derivatives = Path("/Users/dumbeldore/Desktop/derivatives")
    tracts = GenerateTractsDipy(
        derivatives, subj="sub-01", reconstruction="deterministic"
    )
    tracts.run()
