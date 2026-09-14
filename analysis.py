import pandas as pd
import numpy as np
from typing import Tuple, Dict, Any, List, Optional, Sequence
from scipy import stats

from configuration_new import Configuration
from vars import PROFILE_BANDS

default = Configuration()

def cohort_terminal_profile(block_frames: Sequence[pd.DataFrame],
                            cfg: Configuration = default) -> Optional[np.ndarray]:
    """B3: perfil espectral de referencia EXTERNO, mediana da coorte nos blocos finais.

    Calculado uma unica vez, sobre todos os sujeitos, sem qualquer acesso ao rotulo.
    E etapa de COORTE (nao pre-processamento individual) e deve ser descrita como tal.
    """
    rows = []
    for prof in block_frames:
        if len(prof) < 4:
            continue
        k = max(2, min(int(cfg.n_terminal_blocks), len(prof) // 2))
        rows.append(np.median(profile_matrix(prof.iloc[-k:]), axis=0))
    if not rows:
        return None
    ref = np.median(np.vstack(rows), axis=0)
    total = ref.sum()
    return ref / total if total > 0 else ref


def profile_matrix(df: pd.DataFrame) -> np.ndarray:
    return df[["prof_%s" % nm for nm, _, _ in PROFILE_BANDS]].to_numpy(float)

def subject_settling_time(prof: pd.DataFrame, cfg: Configuration = default, cohort_ref_profile: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """B3: instante a partir do qual o registro do sujeito fica ESTAVEL.

    PROBLEMA CORRIGIDO. Na v23 a referencia de estabilidade era a mediana da
    SEGUNDA METADE DO PROPRIO REGISTRO. Sob deriva monotona - e 31 dos 55
    registros tem Spearman(RMS, tempo) < -0,5 - a referencia tambem deriva: o
    metodo compara o sinal com um alvo movel e nao pode, por construcao, detectar
    deriva que atravessa o registro inteiro. Isso tende a declarar estabilidade
    cedo demais.

    CORREÇÃO. A referencia ESPECTRAL passa a ser EXTERNA ao sujeito: mediana da
    coorte nos blocos terminais (``cohort_ref_profile``), calculada uma unica vez
    e cega ao rotulo. O perfil espectral e adimensional (fracoes de banda), logo
    comparavel entre sujeitos.

    LIMITACAO DECLARADA. A referencia de RMS permanece intra-sujeito (mediana dos
    ``n_terminal_blocks`` finais), porque a amplitude em contas depende do ganho
    de contato e NAO e comparavel entre sujeitos. Deriva que abrange todo o
    registro nao e corrigida: e SINALIZADA por ``drift_spans_record``, derivado da
    tendencia de RMS nos blocos terminais. Quando esse sinalizador esta ativo, o
    tempo de acomodacao do sujeito deve ser lido como limite inferior.
    """
    n = len(prof)
    if n < 4:
        return {"settling_time_s": float("nan"), "n_blocks": int(n),
                "settling_status": "insufficient_blocks", "drift_spans_record": False,
                "ref_mode": cfg.settling_ref_mode}
    k = max(2, min(int(cfg.n_terminal_blocks), n // 2))
    tail = prof.iloc[-k:]
    rms_ref = float(np.median(tail["rms_counts"])) + 1e-30
    M = profile_matrix(prof)
    if cfg.settling_ref_mode == "cohort" and cohort_ref_profile is not None:
        prof_ref = np.asarray(cohort_ref_profile, float)
        ref_used = "cohort_terminal_blocks"
    else:
        prof_ref = np.median(profile_matrix(tail), axis=0)
        ref_used = "self_terminal_blocks"
    dev_log2 = np.abs(np.log2((prof["rms_counts"].to_numpy(float) + 1e-30) / rms_ref))
    dev_js = np.array([jensen_shannon(M[i], prof_ref) for i in range(n)])
    stable = (dev_log2 <= cfg.settle_tol_log2) & (dev_js <= cfg.settle_tol_js)
    idx = len(stable)
    for i in range(len(stable) - 1, -1, -1):
        if stable[i]:
            idx = i
        else:
            break
    t = float(prof["t_start_s"].iloc[idx]) if idx < len(prof) else float("nan")
    tail_rho = float(stats.spearmanr(tail["t_start_s"], tail["rms_counts"]).statistic) \
        if k >= 3 else float("nan")
    spans = bool(np.isfinite(tail_rho) and abs(tail_rho) >= cfg.tail_trend_warn)
    return {"settling_time_s": t, "n_blocks": int(n),
            "frac_blocks_stable": float(np.mean(stable)),
            "rms_trend_spearman": float(stats.spearmanr(prof["t_start_s"],
                                                        prof["rms_counts"]).statistic),
            "tail_rms_trend_spearman": tail_rho,
            "drift_spans_record": spans,
            "ref_mode": cfg.settling_ref_mode, "reference_used": ref_used,
            "n_terminal_blocks": int(k),
            "dev_log2_first_block": float(dev_log2[0]), "dev_js_first_block": float(dev_js[0]),
            "settling_status": ("ok" if np.isfinite(t) else "never_stable")}


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Divergência de Jensen-Shannon (base 2) entre dois perfis espectrais em [0,1]."""
    p = np.clip(np.asarray(p, float), 1e-12, None); 
    p = p / p.sum()
    q = np.clip(np.asarray(q, float), 1e-12, None); 
    q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: float(np.sum(a * np.log2(a / b)))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)



def cohort_qc_decision(qc_df: pd.DataFrame, qc: Configuration) -> pd.DataFrame:
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
