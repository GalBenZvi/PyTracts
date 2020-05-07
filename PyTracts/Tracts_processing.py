# import weighted_tracts
import os
from dipy.io.streamline import load_tractogram
import glob
from logs import messages
from pathlib import Path
import time
import nibabel as nib
import logging
from PyTracts import weighted_tracts
from PyTracts.MrtrixTractography import (
    CalculateFibreOrientation,
    GenerateTck,
    GenerateFA,
    GenerateResponses,
    ConvertTck2Trk,
)
from PyTracts.utils import check_dir_existence, check_dir_existence, FSLOUTTYPE
from atlases.atlases import Atlases

ATLAS = Atlases.megaatlas_dir.value


def init_process(mother_dir="/home/gal/Brain_Networks"):
    prep = Prep.Preprocess(mother_dir)
    subjects = prep.subjects
    return mother_dir, subjects


class Generate_Tracts_with_dipy:
    def __init__(self, subject=None):
        self.mother_dir, self.subjects = init_process()
        if subject:
            self.subjects = subject

    def Gen_init_tracts_dipy(self, subj):
        folder_name = f"{self.mother_dir}/Niftis/{subj}"
        (
            gtab,
            data,
            affine,
            labels,
            white_matter,
            nii_file,
            bvec_file,
        ) = weighted_tracts.load_dwi_files(folder_name=folder_name)
        seeds = weighted_tracts.create_seeds_new(labels=labels, affine=affine)
        csd_fit = weighted_tracts.create_csd_model(
            data=data, gtab=gtab, white_matter=white_matter
        )
        fa, classifier = weighted_tracts.create_fa_classifier(
            gtab=gtab, data=data, white_matter=white_matter
        )
        streamlines = weighted_tracts.create_streamlines(
            csd_fit=csd_fit, classifier=classifier, seeds=seeds, affine=affine
        )
        return streamlines, nii_file

    def save_tracts_file(self, subj, streamlines, nii_file):
        folder_name = f"{self.mother_dir}/Derivatives/Streamlines/{subj}"
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
        weighted_tracts.save_ft(folder_name, subj, streamlines, nii_file)

    def run_whole_head_tractography(self):
        for subj in self.subjects:
            streamlines, nii_file = self.Gen_init_tracts_dipy(subj)
            self.save_tracts_file(subj, streamlines, nii_file)


class Generate_Connectivity:
    def __init__(self, derivatives_dir: Path, subj=None, atlas_dir: Path = ATLAS):
        self.derivatives = derivatives_dir
        if subj:
            subjects = [subj]
        else:
            subjects = [subj.name for subj in derivatives_dir.glob("sub-*")]
        subjects.sort()
        subjects_dict = dict()
        for subj in subjects:
            subjects_dict[subj] = derivatives_dir / subj
        self.subjects = subjects_dict
        self.atlas_dir = atlas_dir

    def init_subject_params(self, subj: str):
        dwi_file = (
            self.derivatives
            / subj
            / "dwi"
            / f"{subj}_acq-AP_dwi_preprocessed_biascorr{FSLOUTTYPE}"
        )
        bvec_file = glob.glob(f"{dwi_file.parent}/*.bvec")[0]
        reg_folder = self.derivatives / subj / "atlases"
        stream_folder = self.derivatives / subj / "tractography"
        streamlines_file = stream_folder / "tractogram.trk"
        # streamlines = load_tractogram(streamlines_file, dwi_file)
        streamlines = nib.streamlines.load(streamlines_file)
        return streamlines, stream_folder, reg_folder, bvec_file

    def run_whole_head_connectivity(self):
        for subj in self.subjects:
            (
                streamlines,
                stream_folder,
                reg_folder,
                bvec_file,
            ) = self.init_subject_params(subj)
            lab_labels_index, affine = weighted_tracts.nodes_by_index(reg_folder)
            index_file = glob.glob(f"{self.atlas_dir}/*.txt")[0]
            labels_headers, idx = weighted_tracts.nodes_labels_mega(index_file)
            new_data, m, grouping = weighted_tracts.non_weighted_con_mat_mega(
                streamlines.streamlines, lab_labels_index, affine, idx, stream_folder
            )
            non_weighted_fig_name = (
                f"{stream_folder}/Whole_Head_non-weighted_Connectivity.jpg"
            )
            weighted_tracts.draw_con_mat(
                new_data, labels_headers, non_weighted_fig_name
            )
            weight_by = "1.5_2_AxPasi5"
            weighted_fig_name = (
                f"{stream_folder}/Whole_Head_-{weight_by}_weighted_Connectivity.jpg"
            )
            new_data, mm_weighted = weighted_tracts.weighted_con_mat_mega(
                bvec_file, weight_by, grouping, idx, stream_folder
            )
            weighted_tracts.draw_con_mat(
                new_data, labels_headers, weighted_fig_name, is_weighted=True
            )


class GenerateTractsMrtrix3:
    """
    Mrtrix3-based tractography pipeline
    """

    def __init__(self, mother_dir: Path, subj: str = None):
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

    def __str__(self):
        str_to_print = messages.TRACTS_WITH_MRTRIX3.format(
            mother_dir=self.mother_dir, subjects=self.subjects.keys()
        )
        return str_to_print

    def load_files(self, subj: str, folder_name: Path):
        dwi = Path(
            folder_name / "dwi" / "Mrtrix_prep" / "dwi_preprocessed_biascorr.mif"
        )
        dwi_mask = Path(
            folder_name / "dwi" / "Mrtrix_prep" / "fieldmap_magnitude_brain_mask.mif"
        )
        five_tissue = Path(folder_name / "dwi" / "Mrtrix_prep" / "5TT.mif")
        dwi_nii = Path(
            folder_name / "dwi" / f"{subj}_acq-AP_dwi_preprocessed_biascorr.nii.gz"
        )
        return dwi, dwi_mask, five_tissue, dwi_nii

    def generate_fa(self, dwi: Path, mask: Path, tracts_dir: Path):
        fa = GenerateFA(dwi, mask, tracts_dir)
        return fa

    def generate_responses(self, dwi: Path, mask: Path, tracts_dir: Path):
        resp = GenerateResponses(dwi, mask, tracts_dir)
        return resp

    def calculate_fibre_oriention(
        self, dwi: Path, mask: Path, response_dict: dict, tracts_dir: Path
    ):
        fod = CalculateFibreOrientation(dwi, mask, response_dict, tracts_dir)
        return fod

    def generate_tracts(self, fod_dict: dict, five_tissue: Path, tracts_dir: Path):
        tracts = GenerateTck(fod_dict, five_tissue, tracts_dir)
        return tracts

    def print_start(self, subj: str):
        folder_name = self.subjects[subj]
        tracts_dir = folder_name / "tractography"
        if not tracts_dir.exists():
            tracts_dir.mkdir()
        dwi, dwi_mask, five_tissue, dwi_nii = self.load_files(subj, folder_name)
        str_to_print = messages.MRTRIX_PRINT_START.format(
            subj=subj,
            folder_name=folder_name,
            dwi_name=dwi.name,
            dwi_mask=dwi_mask.name,
            five_tissue=five_tissue.name,
            tracts_dir=tracts_dir.name,
        )
        return tracts_dir, dwi, dwi_mask, five_tissue, dwi_nii, str_to_print

    def convert_tck_to_trk(self, tck_tracts: Path, dwi_nii: Path, trk_tracts: Path):
        trk_converter = ConvertTck2Trk(tck_tracts, dwi_nii, trk_tracts)
        return trk_converter

    def run(self):
        for subj in self.subjects.keys():
            t = time.time()
            (
                tracts_dir,
                dwi,
                mask,
                five_tissue,
                dwi_nii,
                str_to_print,
            ) = self.print_start(subj)
            logging.basicConfig(
                filename=tracts_dir / "tractography.log",
                filemode="w",
                format="%(asctime)s - %(message)s",
                level=logging.INFO,
            )
            logging.info(str_to_print)
            fa = self.generate_fa(dwi, mask, tracts_dir)
            logging.info(fa)
            fa.run()
            resp = self.generate_responses(dwi, mask, tracts_dir)
            logging.info(resp)
            response_dict = resp.run()
            fod = self.calculate_fibre_oriention(dwi, mask, response_dict, tracts_dir)
            logging.info(fod)
            fod_dict = fod.run()
            tract = self.generate_tracts(fod_dict, five_tissue, tracts_dir)
            logging.info(tract)
            tractogram = tract.run()
            trk_file = Path(tractogram.parent / f"{tractogram.stem}.trk")
            trk_converter = self.convert_tck_to_trk(tractogram, dwi_nii, trk_file)
            logging.info(trk_converter)
            tractogram_trk = trk_converter.run()
            elapsed = (time.time() - t) / 60
            logging.info(
                "%s`s whole-brain tractography took %.2f minutes" % (subj, elapsed)
            )


if __name__ == "__main__":
    derivatives = Path("/home/gal/derivatives")
    tracts = GenerateTractsMrtrix3(derivatives, subj="sub-09")
    tracts.run()
