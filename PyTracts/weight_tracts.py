from pathlib import Path
import glob
from PyTract.weighting import weighted_tracts
import nibabel as nib
from PyTracts.utils import FSLOUTTYPE


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
