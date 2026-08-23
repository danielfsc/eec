
@dataclass(frozen=True)
class QCConfig:
    """QC em tres niveis: epoca, participante e coorte. Sempre cego ao rotulo."""
    p2p_ratio_max: float = 5.0
    p2p_ratio_min: float = 0.2
    jump_ratio_max: float = 8.0
    lf_share_ratio_max: float = 2.5
    max_rejected_fraction: float = 0.5
    max_line_noise_ratio: float = 0.25
    max_ocular_index: float = 0.80
    max_muscle_ratio: float = 0.60
    max_zero_diff_fraction: float = 0.20
    max_saturation_fraction: float = 1e-4
    max_nonstationarity_cv: float = 1.00
    min_channel_corr: float = -0.50
    cohort_mad_z_max: float = 4.0
    cohort_outlier_metrics: Tuple[str, ...] = (
        "line_noise_ratio_post", "ocular_index", "muscle_ratio",
        "nonstationarity_cv", "log_rms_counts")
    wrap_fraction_warn: float = 0.60

