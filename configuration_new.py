from dataclasses import dataclass, field
from typing import  Tuple, Optional
import numpy as np
from old_stuff.qc_freeze import QCFreezeConfig
from vars import EXPECTED_SUBJECT_IDS_55

@dataclass(frozen=True)
class Configuration:

    policy_mode: str = "integration"                 # derived | fixed
    ## Acquisition Parameters
    
    bits_resolution = 24

    container_bits: int = 32

    electrode_labels = ['Fp1', 'Fpz', 'Fp2']

    exclusion_time_window = 2.0

    max_data_duration = 7200 # Duração do experimento em segundos

    number_electrodes = 3

    sampling_frequency = 250 #Frequência da coleta de dados em Hz

    fs_verification_tol_hz: float = 0.50

    fs_verification_line_hz: Tuple[float, ...] = (50.0, 60.0)

    fs_is_assumption: bool = True   # B2: revertido para suposicao declarada

    fs_verification_min_prominence_log10: float = 0.30


    # BLOCK BRAKE PARAMETERS
    block_time = 20.0 # em segundos


    #SCHEMA PARAMETERS
    allowed_tasks: Tuple[str, ...] = ("still", "unspecified")

    enforce_allowed_tasks: bool = True

    expected_n_channels: int = 3

    require_integer: bool = True

    require_finite: bool = True


    # #FILTERS PARAMETERS    
    low_filter_frequency = 40.

    high_filter_frequency=1.

    band_limits = [high_filter_frequency,low_filter_frequency]

    filters_order = 4

    notch_quality_factor = 30 #Fator de Qualidade do Filtro Notch em dB

    notch_remove_frequency: Optional[float] = 50.0

    notch_mode: str = "auto"

    notch_line_ratio_threshold: float = 0.02

    # FEATURES PARAMETERS
    ocular_band: Tuple[float, float] = (0.5, 3.0)

    muscle_band: Tuple[float, float] = (20.0, 40.0)

    total_band: Tuple[float, float] = (1.0, 40.0)
    
    bands: Tuple[Tuple[str, float, float], ...] = (
            ("delta", 1.0, 4.0), ("theta", 4.0, 8.0), ("alpha", 8.0, 13.0), ("beta", 13.0, 30.0))

    aperiodic_fit_band: Tuple[float, float] = (3.0, 35.0)

    aperiodic_peak_sd: float = 1.5

    aperiodic_max_iter: int = 5

    alpha_peak_band: Tuple[float, float] = (7.0, 13.0)

    alpha_min_peak_log10: float = 0.05

    perm_entropy_order: int = 3

    perm_entropy_delay: int = 1

    higuchi_kmax: int = 10

    aggregate: str = "median"

    #Temporal PROTOCOL
    fixed_skip_seconds: float = 30.0

    mode: str = "derived"                 # derived | fixed

    fixed_window_seconds: float = 240.0

    min_duration_seconds: float = 600.0

    min_skip_seconds: float = 30.0
    
    max_skip_seconds: float = 420.0

    max_window_seconds: float = 480.0

    n_terminal_blocks: int = 5

    n_sensitivity_windows: int = 3

    settle_quantile: float = 0.80

    settling_ref_mode: str = "cohort"     # B3: cohort | self

    settle_tol_log2: float = 0.35

    settle_tol_js: float = 0.15

    max_within_window_drift_log2: float = 1.50
    
    max_within_window_js: float = 0.35

    tail_trend_warn: float = 0.50
    
    target_window_seconds: float = 240.0

    #PRE PROCESSING PARAMETERS
    epoch_seconds: float = 4.0

    epoch_overlap: float = 0.0

    welch_nperseg_seconds: float = 2.0

    welch_overlap: float = 0.5

    qc_freeze: QCFreezeConfig = field(default_factory=QCFreezeConfig)

    # QC PARAMETERS
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
    # max_nonstationarity_cv: float = 1.00
    max_nonstationarity_cv: float = float("inf")  # B4: desativado
    min_channel_corr: float = -0.50
    cohort_mad_z_max: float = 4.0
    # cohort_outlier_metrics: Tuple[str, ...] = (
    #     "line_noise_ratio_post", "ocular_index", "muscle_ratio",
    #     "nonstationarity_cv", "log_rms_counts")

    cohort_outlier_metrics=(
        "line_noise_ratio_post",
        "ocular_index",
        "muscle_ratio",
        "log_rms_counts",
    ),
    wrap_fraction_warn: float = 0.60#
    min_good_epochs: int = 30



    thresholds_frozen=True,
    freeze_date="2026-10-02",
    freeze_source=(
        "qc_distribution.csv audit: nonstationarity_cv delegado ao drift_pass "
        "para evitar dupla penalizacao na janela de analise. Amostra recuperada: N=32."
    ),
    require_freeze_for_research=True,


    #COHORT PARAMETERS
    expected_n: Optional[int] = 55
    expected_subject_ids: Optional[Tuple[str, ...]] = EXPECTED_SUBJECT_IDS_55
    duplicate_policy: str = "prefer_task"
    preferred_task: str = "still"
    require_one_file_per_subject: bool = True

    #EVIDENCE POLICY PARAMETERS
    min_subjects_for_ml: int = 20
    min_per_class_for_ml: int = 20
    require_both_classes: bool = True
    forbid_clinical_scales: bool = True

    #ANALYSIS PARAMETERS
    covariates: Tuple[str, ...] = ("age", "sex", "education_years")
    fdr_alpha: float = 0.05
    n_boot_effect: int = 2000
    outer_folds: int = 5
    outer_repeats: int = 20
    inner_folds: int = 5
    logistic_c_grid: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0)
    inner_scoring: str = "roc_auc"
    n_boot_metrics: int = 2000
    n_permutations: int = 200
    permutation_outer_repeats: int = 1
    batch_collinearity_warn: float = 0.60
    decision_threshold: float = 0.50            # B11
    permutation_matched_repeats: int = 5        # B1: identico para observado e nulo
    sensitivity_outer_repeats: Optional[int] = None  # B8: None => outer_repeats
    sign_agreement_n_boot: int = 2000
    random_seed: int = 20240517