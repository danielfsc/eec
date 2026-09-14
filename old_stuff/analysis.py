from dataclasses import dataclass
from typing import Tuple, Optional
@dataclass(frozen=True)
class AnalysisConfig:
    """Plano analitico congelado ANTES de ver resultados."""
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

