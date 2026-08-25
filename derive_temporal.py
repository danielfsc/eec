import numpy as np
import pandas as pd
import warnings
from typing import Dict, Any
from temporal_protocol import TemporalProtocolConfig



def derive_temporal_protocol(settle_df: pd.DataFrame, durations: Dict[str, float],
                             tcfg: TemporalProtocolConfig) -> Dict[str, Any]:
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
    skip_raw = float(np.ceil(q / tcfg.block_seconds) * tcfg.block_seconds)
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