import pandas as pd
import numpy as np
import warnings
from typing import Tuple, Dict, Any, List, Optional, Sequence
from scipy import stats

from configuration_new import Configuration
from vars import PROFILE_BANDS, QC_THRESHOLD_MAP

default = Configuration()

# ==============================================================================
# AUDITORIA DE REDUNDÂNCIA: nonstationarity_cv vs drift_pass
# ==============================================================================
def audit_nonstationarity_redundancy(merged: pd.DataFrame, qc_df: pd.DataFrame) -> dict:
    """Investiga se sujeitos estáveis na janela estão sendo eliminados
    apenas pela não-estacionaridade global do registro inteiro."""
    # 1. Recupera as razões individuais de falha no QC
    qc_map = qc_df.set_index("subject_id")["qc_reasons"].to_dict()
    merged["qc_reasons"] = merged["subject_id"].map(qc_map).fillna("")

    # 2. Identifica os 27 sujeitos que passam no drift mas falham no QC
    target_group = merged.loc[(~merged["qc_pass"]) & (merged["drift_pass"])].copy()

    # 3. Categoriza os motivos de falha desse grupo específico
    only_nonstat = []
    nonstat_and_others = []
    other_reasons_only = []

    for _, row in target_group.iterrows():
        sid = row["subject_id"]
        reasons = [r.strip() for r in row["qc_reasons"].split(";") if r.strip()]
        has_nonstat = any(r.startswith("nonstationarity") or "z__nonstationarity_cv" in r for r in reasons)
        other_reasons = [r for r in reasons if not (r.startswith("nonstationarity") or "z__nonstationarity_cv" in r)]
        
        if has_nonstat and len(other_reasons) == 0:
            only_nonstat.append(sid)
        elif has_nonstat and len(other_reasons) > 0:
            nonstat_and_others.append((sid, other_reasons))
        else:
            other_reasons_only.append((sid, other_reasons))

    redundancy_report = {
        "n_target_qc_fail_drift_pass": len(target_group),
        "n_failed_ONLY_by_nonstationarity": len(only_nonstat),
        "ids_failed_ONLY_by_nonstationarity": only_nonstat,
        "n_failed_by_nonstat_plus_others": len(nonstat_and_others),
        "n_failed_by_other_qc_reasons": len(other_reasons_only),
        "potential_analytic_n": int(merged["qc_pass"].sum() + len(only_nonstat)),
    }
    
    print("\n===== AUDITORIA DE REDUNDÂNCIA (QC vs DRIFT) =====")
    print(f"Sujeitos com drift_pass=True mas qc_pass=False: {len(target_group)}")
    print(f"  -> Falham ESTRITAMENTE por nonstationarity_cv: {len(only_nonstat)}")
    print(f"  -> Falham por nonstationarity + outros motivos: {len(nonstat_and_others)}")
    print(f"  -> Falham apenas por outros motivos de QC:     {len(other_reasons_only)}")
    print(f"N amostral se nonstationarity for delegado ao drift: {len(merged.loc[merged['drift_pass']]) - len(other_reasons_only)}")
    print("==================================================\n")
    
    return redundancy_report

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

def derive_temporal_protocol(settle_df: pd.DataFrame, durations: Dict[str, float],
                             tcfg: Configuration) -> Dict[str, Any]:
    """A5/B4: converte os tempos de acomodacao individuais em UMA regra de coorte.

    Passos: (i) elegibilidade por duracao minima pre-declarada; (ii) inicio comum =
    quantil dos tempos de acomodacao, na grade de blocos, limitado por
    [min_skip_seconds, max_skip_seconds]; (iii) duracao comum = menor sobra entre os
    elegiveis, limitada por target/max; (iv) janelas de sensibilidade nao sobrepostas.
    Nenhum rotulo participa de qualquer etapa.

    B4 CORRIGIDO. Na v23, se o quantil dos tempos de acomodacao excedesse
    ``max_skip_seconds``, o inicio era truncado EM SILENCIO e a janela podia voltar
    a cair no transitorio - exatamente o problema que a derivacao existe para
    evitar. Agora o truncamento e explicito (``skip_capped``, ``window_capped``),
    emite ``warnings.warn`` e marca ``protocol_valid=False``, para que o resultado
    nao seja lido como se a regra tivesse sido respeitada.
    """
    if tcfg.mode == "fixed":
        skip, win = tcfg.fixed_skip_seconds, tcfg.fixed_window_seconds
        eligible = [s for s, d in durations.items() if d >= skip + win]
        return {"rule": "fixed", "skip_seconds": skip, "window_seconds": win,
                "eligible_subjects": sorted(eligible),
                "excluded_short": sorted(set(durations) - set(eligible)),
                "sensitivity_windows": [(skip, win)], "settling_quantile_used": None,
                "skip_capped": False, "window_capped": False, "protocol_valid": True,
                "n_drift_spans_record": 0}
    eligible = [s for s, d in durations.items() if d >= tcfg.min_duration_seconds]
    excluded = sorted(set(durations) - set(eligible))
    if not eligible:
        raise ValueError("Nenhum sujeito atinge min_duration_seconds=%.0f s."
                         % tcfg.min_duration_seconds)
    sdf = settle_df.set_index("subject_id").reindex(eligible)
    st = pd.to_numeric(sdf["settling_time_s"], errors="coerce").to_numpy(float)
    st = st[np.isfinite(st)]
    n_spans = int(pd.Series(sdf.get("drift_spans_record", pd.Series(dtype=bool))).fillna(False).sum())
    q = float(np.quantile(st, tcfg.settle_quantile)) if st.size else tcfg.min_skip_seconds
    skip_raw = float(np.ceil(q / tcfg.block_time) * tcfg.block_time)
    skip = float(min(max(skip_raw, tcfg.min_skip_seconds), tcfg.max_skip_seconds))
    skip_capped = bool(skip_raw > tcfg.max_skip_seconds)
    if skip_capped:
        warnings.warn("B4: inicio derivado (%.0f s) excede max_skip_seconds (%.0f s) e foi "
                      "truncado; a janela pode conter transitorio de acomodacao."
                      % (skip_raw, tcfg.max_skip_seconds), RuntimeWarning)
    avail = min(durations[s] for s in eligible) - skip
    win_raw = min(avail, tcfg.target_window_seconds)
    win = float(np.floor(min(win_raw, tcfg.max_window_seconds)))
    window_capped = bool(avail < tcfg.target_window_seconds)
    if window_capped:
        warnings.warn("B4: duracao comum limitada pela sobra disponivel (%.0f s) e nao pelo "
                      "alvo (%.0f s)." % (avail, tcfg.target_window_seconds), RuntimeWarning)
    if win <= 0:
        raise ValueError("Janela comum nao positiva; revise min_duration_seconds.")
    k = max(1, int(np.floor(avail / win)))
    k = min(k, tcfg.n_sensitivity_windows)
    if k > 1:
        span = avail - win
        starts = [float(np.floor(skip + span * i / (k - 1))) for i in range(k)]
        for i in range(1, k):  # garante nao sobreposicao
            starts[i] = max(starts[i], starts[i - 1] + win)
        starts = [s for s in starts if s + win <= skip + avail + 1e-9]
    else:
        starts = [skip]
    return {"rule": "derived_from_data", "skip_seconds": skip, "window_seconds": win,
            "skip_seconds_before_cap": skip_raw,
            "skip_capped": skip_capped, "window_capped": window_capped,
            "protocol_valid": bool(not skip_capped),
            "settling_quantile_used": tcfg.settle_quantile,
            "settling_time_quantile_s": q,
            "cohort_settling_median_s": float(np.median(st)) if st.size else None,
            "cohort_settling_max_s": float(np.max(st)) if st.size else None,
            "available_after_skip_s": float(avail),
            "n_drift_spans_record": n_spans,
            "settling_reference_mode": tcfg.settling_ref_mode,
            "eligible_subjects": sorted(eligible), "excluded_short": excluded,
            "sensitivity_windows": [(float(s), float(win)) for s in starts]}

def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Divergência de Jensen-Shannon (base 2) entre dois perfis espectrais em [0,1]."""
    #Garantindo que p não seja 0;
    p = np.clip(np.asarray(p, float), 1e-12, None); 
    p = p / p.sum()
    q = np.clip(np.asarray(q, float), 1e-12, None); 
    q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: float(np.sum(a * np.log2(a / b)))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

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
    # Sem blocos o suficiente
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
    #Por que log2?
    dev_log2 = np.abs(np.log2((prof["rms_counts"].to_numpy(float) + 1e-30) / rms_ref))
    dev_js = np.array([jensen_shannon(M[i], prof_ref) for i in range(n)])
    stable = (dev_log2 <= cfg.settle_tol_log2) & (dev_js <= cfg.settle_tol_js)
    # Aqui está com um problema: Alguns sujeitos tem INSTABILIDADE em, por exemplo, idx=58! Problema na divergência do Log2
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

def qc_distribution_report(qc_df: pd.DataFrame, qc: Configuration) -> pd.DataFrame:
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
    # print(qc.cohort_outlier_metrics)
    for m in qc.cohort_outlier_metrics:
        # print(f"Processing cohort outlier metric: {m}")
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




