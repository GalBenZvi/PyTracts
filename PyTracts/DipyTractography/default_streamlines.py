from dipy.data import default_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from pathlib import Path
import numpy as np


class DefaultStreamlines:
    def __init__(
        self,
        folder_name: Path,
        csa_model,
        data: np.ndarray,
        white_matter: np.ndarray,
        relative_peak_threshold: float = 0.8,
        min_angle: float = 45,
    ):
        self.folder_name = folder_name
        self.csa_model = csa_model
        self.data = data
        self.white_matter = white_matter
        self.relative_peak_threshold = relative_peak_threshold
        self.min_angle = min_angle

    def generate_peaks(self):
        csa_peaks = peaks_from_model(
            self.csa_model,
            self.data,
            default_sphere,
            relative_peak_threshold=self.relative_peak_threshold,
            min_separation_angle=self.min_angle,
            mask=self.white_matter,
        )
        return csa_peaks
