from dataclasses import dataclass
from typing import Tuple, List
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

import pandas as pd
import numpy as np
from vars import QC_THRESHOLD_MAP
from qc_freeze import QCFreezeConfig
from evidence_policy import EvidencePolicy

def assert_qc_thresholds_frozen(freeze: QCFreezeConfig, policy: EvidencePolicy) -> None:
    """B5: trava dura - modo research exige limiares de QC congelados e datados."""
    if policy.mode != "research" or not freeze.require_freeze_for_research:
        return
    if not freeze.thresholds_frozen:
        raise RuntimeError(
            "B5: modo research bloqueado. Os limiares de QC ainda nao foram congelados. "
            "Execute em modo integration, inspecione 02_quality/qc_distribution.csv (que "
            "informa quantos sujeitos cada limiar excluiria), fixe os valores e declare "
            "QCFreezeConfig(thresholds_frozen=True, freeze_date=..., freeze_source=...).")


def cohort_qc_decision(qc_df: pd.DataFrame, qc: QCConfig) -> pd.DataFrame:
    """Decisao de exclusao PRE-ESPECIFICADA: limiares absolutos + outlier de coorte (MAD)."""
    out = qc_df.copy()
    reasons: List[List[str]] = [[] for _ in range(len(out))]
    def _abs_rule(col: str, thr: float, greater: bool, tag: str) -> None:
        if col not in out.columns:
            return
        v = pd.to_numeric(out[col], errors="coerce").to_numpy(float)
        bad = (v > thr) if greater else (v < thr)
        for i in np.where(np.isfinite(v) & bad)[0]:
            reasons[i].append("%s=%.4g" % (tag, v[i]))

    _abs_rule("line_noise_ratio_post", qc.max_line_noise_ratio, True, "line_noise")
    _abs_rule("ocular_index", qc.max_ocular_index, True, "ocular")
    _abs_rule("muscle_ratio", qc.max_muscle_ratio, True, "muscle")
    _abs_rule("zero_diff_fraction", qc.max_zero_diff_fraction, True, "frozen_samples")
    _abs_rule("saturation_fraction", qc.max_saturation_fraction, True, "saturation")
    _abs_rule("nonstationarity_cv", qc.max_nonstationarity_cv, True, "nonstationarity")
    _abs_rule("min_channel_corr", qc.min_channel_corr, False, "channel_corr")
    if "flat_channels" in out.columns:
        for i in np.where(pd.to_numeric(out["flat_channels"], errors="coerce").fillna(0) > 0)[0]:
            reasons[i].append("flat_channel")
    # print(out)
    for m in qc.cohort_outlier_metrics:
        if m not in out.columns:
            continue
        v = pd.to_numeric(out[m], errors="coerce").to_numpy(float)
        med = np.nanmedian(v)
        mad = np.nanmedian(np.abs(v - med))
        scale = 1.4826 * mad if mad > 0 else np.nanstd(v)
        z = (v - med) / (scale + 1e-30)
        out["z__" + m] = z
        for i in np.where(np.isfinite(z) & (np.abs(z) > qc.cohort_mad_z_max))[0]:
            reasons[i].append("cohort_outlier_%s(z=%.1f)" % (m, z[i]))
    out["qc_reasons"] = ["; ".join(r) for r in reasons]
    out["qc_pass"] = [len(r) == 0 for r in reasons]
    return out

def qc_distribution_report(qc_df: pd.DataFrame, qc: QCConfig) -> pd.DataFrame:
    """B5: confronta CADA limiar de QC com a distribuicao empirica da coorte.

    ERRO CORRIGIDO. Na v23 os limiares (max_ocular_index=0,80, max_line_noise_ratio
    =0,25, max_muscle_ratio=0,60, max_nonstationarity_cv=1,00) nunca foram
    comparados com os valores realmente observados. O numero de exclusoes que eles
    produzem era desconhecido a priori e podia variar de zero a quase toda a coorte
    - e ninguem saberia antes de rodar.

    Esta tabela informa, por metrica: n valido, minimo, quartis, maximo, o limiar
    vigente e QUANTOS sujeitos ele excluiria (contagem e fracao). E o insumo
    obrigatorio para congelar os limiares antes de habilitar o modo research.
    Nao usa rotulo em nenhum momento.
    """
    rows = []
    for col, (attr, direction) in QC_THRESHOLD_MAP.items():
        if col not in qc_df.columns:
            continue
        v = pd.to_numeric(qc_df[col], errors="coerce").to_numpy(float)
        fin = v[np.isfinite(v)]
        thr = float(getattr(qc, attr))
        n_ex = int(np.sum(fin > thr) if direction == "greater" else np.sum(fin < thr))
        rows.append({"metric": col, "threshold_param": attr, "threshold_value": thr,
                     "direction": direction, "n_valid": int(fin.size),
                     "min": float(fin.min()) if fin.size else np.nan,
                     "p05": float(np.percentile(fin, 5)) if fin.size else np.nan,
                     "q25": float(np.percentile(fin, 25)) if fin.size else np.nan,
                     "median": float(np.median(fin)) if fin.size else np.nan,
                     "q75": float(np.percentile(fin, 75)) if fin.size else np.nan,
                     "p95": float(np.percentile(fin, 95)) if fin.size else np.nan,
                     "max": float(fin.max()) if fin.size else np.nan,
                     "n_excluded_by_threshold": n_ex,
                     "frac_excluded_by_threshold": float(n_ex / fin.size) if fin.size else np.nan})
    for m in qc.cohort_outlier_metrics:
        col = "z__" + m
        if col not in qc_df.columns:
            continue
        z = pd.to_numeric(qc_df[col], errors="coerce").to_numpy(float)
        fin = z[np.isfinite(z)]
        n_ex = int(np.sum(np.abs(fin) > qc.cohort_mad_z_max))
        rows.append({"metric": col, "threshold_param": "cohort_mad_z_max",
                     "threshold_value": float(qc.cohort_mad_z_max), "direction": "abs_greater",
                     "n_valid": int(fin.size),
                     "min": float(fin.min()) if fin.size else np.nan,
                     "p05": float(np.percentile(fin, 5)) if fin.size else np.nan,
                     "q25": float(np.percentile(fin, 25)) if fin.size else np.nan,
                     "median": float(np.median(fin)) if fin.size else np.nan,
                     "q75": float(np.percentile(fin, 75)) if fin.size else np.nan,
                     "p95": float(np.percentile(fin, 95)) if fin.size else np.nan,
                     "max": float(fin.max()) if fin.size else np.nan,
                     "n_excluded_by_threshold": n_ex,
                     "frac_excluded_by_threshold": float(n_ex / fin.size) if fin.size else np.nan})
    return pd.DataFrame(rows)
