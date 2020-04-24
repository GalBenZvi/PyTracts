import nipype.interfaces.mrtrix3 as mrt
import os
import nibabel as nib
from pathlib import Path
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes


def fit_tensors(dwi_file: Path, mask_file: Path, dti_file: Path, fa_file: Path):
    dilated_mask = mask_file.parent / f"{mask_file.stem}_dilated.mif"
    cmd = f"maskfilter {mask_file} dilate {dilated_mask} -npass 3"
    os.system(cmd)
    tsr = mrt.FitTensor()
    tsr.inputs.in_file = dwi_file
    tsr.inputs.in_mask = dilated_mask
    tsr.inputs.out_file = dti_file
    print(tsr.cmdline)
    tsr.run()
    comp = mrt.TensorMetrics()
    comp.inputs.in_file = dti_file
    comp.inputs.out_fa = fa_file
    print(comp.cmdline)
    comp.run()
    dti_file.unlink()
    return comp.inputs.out_fa
    # tsr.inputs.args = '- | tensor2metrtic - -fa'


def gen_response(dwi_file: Path, dwi_mask: Path, working_dir: Path):
    resp = mrt.ResponseSD()
    resp.inputs.in_file = dwi_file
    resp.inputs.algorithm = "dhollander"
    resp.inputs.csf_file = working_dir / "response_csf.txt"
    resp.inputs.wm_file = working_dir / "response_wm.txt"
    resp.inputs.gm_file = working_dir / "response_gm.txt"
    resp.inputs.in_mask = dwi_mask
    print(resp.cmdline)
    resp.run()
    return resp.inputs.wm_file, resp.inputs.gm_file, resp.inputs.csf_file


def calc_fibre_orientation(
    dwi_file: str, dwi_mask: str, response_dict: dict, fod_dict: dict
):
    fod = mrt.EstimateFOD()
    fod.inputs.algorithm = "msmt_csd"
    fod.inputs.in_file = dwi_file
    fod.inputs.wm_txt = response_dict["wm"]
    fod.inputs.wm_odf = fod_dict["wm"]
    fod.inputs.gm_txt = response_dict["gm"]
    fod.inputs.gm_odf = fod_dict["gm"]
    fod.inputs.csf_txt = response_dict["csf"]
    fod.inputs.csf_odf = fod_dict["csf"]
    fod.inputs.mask_file = dwi_mask
    fod.inputs.max_sh = 10, 0, 0
    print(fod.cmdline)
    fod.run()


def generate_tracts(fod_wm: Path, tractogram: Path, seg_5tt: Path):
    tk = mrt.Tractography()
    tk.inputs.in_file = fod_wm
    tk.inputs.out_file = tractogram
    tk.inputs.act_file = seg_5tt
    tk.inputs.backtrack = True
    tk.inputs.algorithm = "iFOD2"
    tk.inputs.max_length = 250
    tk.inputs.out_seeds = tractogram.parent / "seeds.csv"
    # tk.inputs.power = 3
    tk.inputs.crop_at_gmwmi = True
    tk.inputs.seed_dynamic = fod_wm
    tk.inputs.select = 350000
    print(tk.cmdline)
    tk.run()
    return tk.inputs.out_file


def convert_tck_to_trk(tracts_tck: Path, dwi_file: Path, tracts_trk: Path):
    nii = nib.load(dwi_file)
    header = {}
    header[Field.VOXEL_TO_RASMM] = nii.affine.copy()
    header[Field.VOXEL_SIZES] = nii.header.get_zooms()[:3]
    header[Field.DIMENSIONS] = nii.shape[:3]
    header[Field.VOXEL_ORDER] = "".join(aff2axcodes(nii.affine))
    tck = nib.streamlines.load(tracts_tck)
    nib.streamlines.save(tck.tractogram, str(tracts_trk), header=header)
    return tracts_trk
