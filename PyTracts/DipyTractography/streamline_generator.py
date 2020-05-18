from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.streamline import Streamlines
from dipy.data import default_sphere, small_sphere
from dipy.direction import ProbabilisticDirectionGetter
from dipy.direction import DeterministicMaximumDirectionGetter
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import save_trk
from logs import messages
import numpy as np
from pathlib import Path
import nibabel as nib


class StreamlineGenerator:
    def __init__(
        self,
        folder_name: Path,
        csd_fit,
        hardi_img: nib.Nifti1Image,
        reconstruction: str,
        stopping_criterion,
        seeds,
        affine: np.ndarray,
        tractogram_fname: str = "Whole_brain_tractography.trk",
        sphere: str = "default",
        max_angle: float = 30.0,
        step_size: float = 0.5,
        csa_peaks=None,
    ):
        self.folder_name = folder_name
        self.reconstruction = reconstruction
        self.stopping_criterion = stopping_criterion
        self.csd_fit = csd_fit
        self.affine = affine
        self.hardi_img = hardi_img
        self.tractogram_fname = tractogram_fname
        if sphere == "small":
            self.sphere = small_sphere
        elif sphere == "default":
            self.sphere = default_sphere
        self.seeds = seeds
        self.step_size = step_size
        self.max_angle = max_angle
        self.csa_peaks = csa_peaks

    def __str__(self):
        str_to_print = messages.GENERATE_STREAMLINES(
            folder=self.folder_name,
            recon=self.reconstruction,
            fname=self.tractogram_fname,
            sphere=self.sphere,
            angle=self.max_angle,
            step=self.step_size,
            tracts_loc=self.folder_name / self.tractogram_fname,
        )
        return str_to_print

    def get_directions(self):
        if self.reconstruction.lower() == "deterministic":
            directions_getter = DeterministicMaximumDirectionGetter.from_shcoeff(
                self.csd_fit.shm_coeff, max_angle=self.max_angle, sphere=self.sphere
            )
        elif self.reconstruction.lower() == "probabilistic":
            directions_getter = ProbabilisticDirectionGetter.from_shcoeff(
                self.csd_fit.shm_coeff, max_angle=self.max_angle, sphere=self.sphere
            )
        elif self.reconstruction.lower() == "default":
            directions_getter = self.csa_peaks
        return directions_getter

    def generate_streamlines(self, directions_getter):
        streamline_generator = LocalTracking(
            directions_getter,
            self.stopping_criterion,
            self.seeds,
            self.affine,
            step_size=self.step_size,
        )
        streamlines = Streamlines(streamline_generator)
        return streamlines

    def save_streamlines(self, streamlines):
        tractogram = StatefulTractogram(streamlines, self.hardi_img, Space.RASMM)
        save_trk(tractogram, str(self.folder_name / self.tractogram_fname))

    def run(self):
        directions_getter = self.get_directions()
        streamlines = self.generate_streamlines(directions_getter)
        self.save_streamlines(streamlines)
        return streamlines
