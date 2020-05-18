from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, auto_response
from dipy.reconst.shm import CsaOdfModel
from dipy.data import default_sphere
from dipy.direction import peaks_from_model
from dipy.core.gradients import GradientTable
import dipy.reconst.dti as dti
from dipy.reconst.dti import fractional_anisotropy
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from logs import messages


class CreateModel:
    def __init__(
        self,
        folder_name: Path,
        data: np.ndarray,
        gtab: GradientTable,
        white_matter: np.ndarray,
        sh_order: int = 6,
        roi_radius: int = 10,
        fa_thr: float = 0.7,
        relative_peak_threshold: float = 0.8,
        min_separation_angle: int = 45,
        stopping_threshold: float = 0.25,
    ):
        self.folder_name = folder_name
        self.stopping_threshold = stopping_threshold
        self.data = data
        self.gtab = gtab
        self.white_mask = white_matter
        self.sh_order = sh_order
        self.roi_radius = roi_radius
        self.fa_thr = fa_thr
        self.relative_peak_threshold = relative_peak_threshold
        self.min_separation_angle = min_separation_angle

    def __str__(self):
        str_to_print = messages.CREATE_MODEL(
            sh=self.sh_order,
            roi=self.roi_radius,
            fa=self.fa_thr,
            relative=self.relative_peak_threshold,
            min_separation=self.min_separation_angle,
        )
        return str_to_print

    def generate_response(self):
        response, ratio = auto_response(
            self.gtab, self.data, roi_radius=self.roi_radius, fa_thr=self.fa_thr
        )
        return response, ratio

    def calculate_model(self, response):
        csa_model = CsaOdfModel(self.gtab, sh_order=self.sh_order)
        csd_model = ConstrainedSphericalDeconvModel(
            self.gtab, response, sh_order=self.sh_order
        )
        csd_fit = csd_model.fit(self.data, mask=self.white_mask)
        return csd_fit, csa_model

    def create_stopping_criterion(self, csa_model):
        tensor_model = dti.TensorModel(self.gtab)
        tenfit = tensor_model.fit(self.data, mask=self.white_mask)
        fa = fractional_anisotropy(tenfit.evals)
        stopping_criterion = ThresholdStoppingCriterion(fa, 0.18)
        # gfa = csa_model.fit(self.data, mask=self.white_mask).gfa
        # stopping_criterion = ThresholdStoppingCriterion(gfa, self.stopping_threshold)
        return stopping_criterion

    def quality_assurance(self, model):
        data = model.gfa
        sli = data.shape[2] // 2
        fig, axs = plt.subplots(1, 2, constrained_layout=True)
        title_obj = fig.suptitle(
            "Generalized fractional anisotrophy (GFA) quality assurance"
        )
        plt.setp(title_obj, color="w")
        titles = ["GFA", "Tracking mask"]
        plot_data = [data[:, :, sli].T, (data[:, :, sli] > 0.25).T]
        for ax, cur_plot, title in zip(axs, plot_data, titles):
            ax.imshow(cur_plot, cmap="gray", origin="lower")
            ax.set_axis_off()
            title_obj = ax.set_title(title)
            plt.setp(title_obj, color="w")
        plt.savefig(
            Path(self.folder_name / "GFA_quality_assurance.png"),
            bbox_inches="tight",
            facecolor="black",
        )

    def run(self):
        response, ratio = self.generate_response()
        odf_model, csa_model = self.calculate_model(response)
        self.quality_assurance(odf_model)
        stopping_criterion = self.create_stopping_criterion(csa_model)
        return odf_model, stopping_criterion, csa_model
