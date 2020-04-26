# import weighted_tracts
import os
from dipy.io.streamline import load_tractogram
import glob
from PyTracts.code import Mrtrix3_methods as mrt_methods
from PyTracts.code import messages
from pathlib import Path
import time


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
    def __init__(self, subject=None):
        self.mother_dir, self.subjects = init_process()
        if subject:
            self.subjects = [subject]

    def init_subject_params(self, subj: str):
        dwi_file = f"{self.mother_dir}/Niftis/{subj}/dwi/diff_corrected.nii.gz"
        bvec_file = glob.glob(f"{self.mother_dir}/Niftis/{subj}/dwi/*.bvec")[0]
        reg_folder = (
            f"{self.mother_dir}/Derivatives/Registrations/{subj}/Atlases_and_Transforms"
        )
        stream_folder = f"{self.mother_dir}/Derivatives/Streamlines/{subj}"
        streamlines_file = f"{stream_folder}/{subj}_wholebrain.trk"
        streamlines = load_tractogram(streamlines_file, dwi_file)

        return streamlines, stream_folder, reg_folder, bvec_file

    def run_whole_head_connectivity(self):
        atlas_folder = f"{self.mother_dir}/Derivatives/megaatlas"
        for subj in self.subjects:
            (
                streamlines,
                stream_folder,
                reg_folder,
                bvec_file,
            ) = self.init_subject_params(subj)
            lab_labels_index, affine = weighted_tracts.nodes_by_index(reg_folder)
            index_file = f"{atlas_folder}/megaatlascortex2nii_origin.txt"
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


def check_files_existence(files: list):
    """
    Checks the existence of a list of files (all should exist)
    Arguments:
        file {Path} -- [description]

    Returns:
        [type] -- [description]
    """
    exist = all(Path(f).exists() for f in files)
    return exist


def check_dir_existence(directory: Path):
    if not directory.exists():
        print(f"Creating directory: {directory}")
        directory.mkdir()


class GenerateFA:
    def __init__(self, dwi_file: Path, mask_file: Path, tracts_dir: Path):
        """
        Class to generate FA and DTI images from preprocessed dwi image and mask.
        Arguments:
            tracts_dir {Path} -- [output directory, where fa.mif and dti.if will be produced]
            dwi_file {Path} -- [path to preprocessed dwi image]
            mask_file {Path} -- [path to dwi mask image]
        """
        self.tract_dir = tracts_dir
        check_dir_existence(self.tract_dir)
        self.dwi_file = dwi_file
        self.mask_file = mask_file
        self.dti_file = tracts_dir / "dti.mif"
        self.fa_file = tracts_dir / "fa.mif"
        self.exist = check_files_existence([self.fa_file, self.dti_file])

    def __str__(self):
        str_to_print = messages.GENERATE_FA.format(
            subj_dir=self.tract_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask_file.name,
            tracts_dir=self.tract_dir,
            FA_name=self.fa_file.name,
        )
        return str_to_print

    def generate_fa(self):
        fa_file = mrt_methods.fit_tensors(
            self.dwi_file, self.mask_file, self.dti_file, self.fa_file
        )
        return fa_file

    def run(self):
        if not self.exist:
            print("Generating FA image for group-level analysis")
            self.fa_file = self.generate_fa()
        else:
            print("Already generated FA image for group-level analysis, continuing...")
        return self.fa_file


class GenerateResponses:
    def __init__(self, dwi_file: Path, mask: Path, tracts_dir: Path):
        """
        Estimating tissue response functions for spherical deconvolution
        Arguments:
            dwi_file {Path} -- [path to preprocessed dwi image]
            mask {Path} -- [path to dwi mask image]
            tracts_dir {Path} -- [path to output directory, where response_*tissue*.txt file will be produced]

        Returns:
            [type] -- [description]
        """
        self.dwi_file = dwi_file
        self.mask = mask
        self.tracts_dir = tracts_dir
        self.response_wm, self.response_gm, self.response_csf = [
            tracts_dir / f
            for f in ["response_wm.txt", "response_gm.txt", "response_csf.txt"]
        ]
        self.exist = check_files_existence(
            [self.response_wm, self.response_gm, self.response_csf]
        )

    def __str__(self):
        str_to_print = messages.GENERATE_RESPONSES.format(
            subj_dir=self.tracts_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask.name,
            tracts_dir=self.tracts_dir,
        )
        return str_to_print

    def generate_response(self):
        response_wm, response_gm, response_csf = mrt_methods.generate_response(
            self.dwi_file, self.mask, self.tracts_dir
        )
        return response_wm, response_gm, response_csf

    def run(self):
        if not self.exist:
            print("Estimating tissue response functions for spherical deconvolution")
            self.response_wm, self.response_gm, self.response_csf = (
                self.generate_response()
            )
        else:
            print("Already estimated tissue response functions, continuing...")
        response_dict = dict()
        for key, val in zip(
            ["wm", "gm", "csf"], [self.response_wm, self.response_gm, self.response_csf]
        ):
            response_dict[key] = val
        return response_dict


class CalculateFibreOrientation:
    def __init__(
        self, dwi_file: Path, dwi_mask: Path, response_dict: dict, tracts_dir: Path
    ):
        """
        Estimation of fibre orientation distributions
        Arguments:
            dwi_file {Path} -- [path to preprocessed dwi file]
            response_dict {dict} -- [dictionary with "wm","gm","csf" as keys, and paths to corresponding response_{tissue}.txt files as values]

        Returns:
            [type] -- [description]
        """
        self.dwi_file = dwi_file
        self.tracts_dir = tracts_dir
        self.mask = dwi_mask
        self.response_dict = response_dict
        self.fod_dict = dict()
        for tissue in ["wm", "gm", "csf"]:
            self.fod_dict[tissue] = tracts_dir / f"FOD_{tissue}.mif"
        self.exist = check_files_existence(list(self.fod_dict.values()))

    def __str__(self):
        str_to_print = messages.CALCULATE_FIBER_ORIENTATION.format(
            subj_dir=self.tracts_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask.name,
            tracts_dir=self.tracts_dir,
        )
        return str_to_print

    def fibre_orientation(self):
        mrt_methods.calculate_fibre_orientation(
            self.dwi_file, self.mask, self.response_dict, self.fod_dict
        )

    def run(self):
        if not self.exist:
            print("Estimating Fibre Orientation Distributions")
            self.fibre_orientation()
        else:
            print("Already estimated fibre orientation distributions, continuing...")
        return self.fod_dict


class GenerateTck:
    def __init__(self, fod_dict: dict, seg_5tt: Path, tract_dir: Path):
        """
        Generate .tck tracts file using Mrtrix3's iFOD2 algorithm
        Arguments:
            fod_dict {dict} -- [dictionary with "wm","gm","csf" as keys, and paths to corresponding FOD_{tissue}.mif files as values]
            tract_dir {Path} -- [Path to directory containing all tracts-processing-related files]
            seg_5tt {Path} -- [Path to 5-tissue-type.mif file]

        Returns:
            [type] -- [description]
        """
        self.fod_wm = fod_dict["wm"]
        self.seg_5tt = seg_5tt
        self.tractogram = tract_dir / "tractogram.tck"
        self.exist = check_files_existence([self.tractogram])

    def __str__(self):
        str_to_print = messages.GENERATE_TCK.format(
            subj_dir=self.tractogram.parent.parent,
            fod_wm_name=self.fod_wm.name,
            tracts_dir=self.tractogram.parent,
            seg_5tt_name=self.seg_5tt.name,
            tractogram_name=self.tractogram.name,
        )
        return str_to_print

    def generate_tracts(self):
        tractogram = mrt_methods.generate_tracts(
            self.fod_wm, self.tractogram, self.seg_5tt
        )
        print(f"Generated tractogram.tck file at {tractogram.parent}")

    def run(self):
        if not self.exist:
            print("Generating tractogram.tck file using iFOD2 algorithm...")
            self.generate_tracts()
        else:
            print(
                "Already generated tractogram.\n to recreate it, please remove the currently existing .tck file"
            )
        return self.tractogram


class ConvertTck2Trk:
    def __init__(self, tck_stream: Path, dwi_nii: Path, trk_stream: Path):
        self.tck_stream = tck_stream
        self.dwi = dwi_nii
        self.trk_stream = trk_stream
        self.exist = check_files_existence([self.trk_stream])

    def __str__(self):
        str_to_print = messages.CONVERT_TCK_2_TRK.format(
            subj_dir=self.dwi.parent.parent,
            dwi_name=self.dwi.name,
            tck_name=self.tck_stream.name,
            trk_name=self.trk_stream.name,
        )
        return str_to_print

    def convert(self):
        trk_file = mrt_methods.convert_tck_to_trk(
            self.tck_stream, self.dwi, self.trk_stream
        )

    def run(self):
        if not self.exist:
            print("Converting tractography file from .tck format to .trk")
            self.convert()
        else:
            print("Given input for .trk tractography file already exists.")
        return self.trk_stream


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
        dwi, dwi_mask, five_tissue, dwi_nii = self.load_files(subj, folder_name)
        str_to_print = messages.MRTRIX_PRINT_START.format(
            subj=subj,
            folder_name=folder_name,
            dwi_name=dwi.name,
            dwi_mask=dwi_mask.name,
            five_tissue=five_tissue.name,
            tract_dir=tracts_dir.name,
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
            print(str_to_print)
            fa = self.generate_fa(dwi, mask, tracts_dir)
            print(fa)
            fa.run()
            resp = self.generate_responses(dwi, mask, tracts_dir)
            print(resp)
            response_dict = resp.run()
            fod = self.calculate_fibre_oriention(dwi, mask, response_dict, tracts_dir)
            print(fod)
            fod_dict = fod.run()
            tract = self.generate_tracts(fod_dict, five_tissue, tracts_dir)
            print(tract)
            tractogram = tract.run()
            trk_file = Path(tractogram.parent / f"{tractogram.stem}.trk")
            trk_converter = self.convert_tck_to_trk(tractogram, dwi_nii, trk_file)
            print(trk_converter)
            tractogram_trk = trk_converter.run()
            elapsed = (time.time() - t) / 60
            print("%s`s whole-brain tractography took %.2f minutes" % (subj, elapsed))
