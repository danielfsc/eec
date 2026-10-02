import warnings
import math
import hashlib
import json
import numpy as np
import pandas as pd
import filter as ft
import scipy.signal as sg
import scipy.stats as stats
import analysis as an
from typing import Tuple, Dict, Any, List, Sequence, Optional
from pathlib import Path
from dataclasses import asdict
from vars import PROFILE_BANDS, NAME_REGEX, PRIMARY_FEATURES, CLINICAL_SCALE_TOKENS
import filter as ft
from configuration_new import Configuration
from old_stuff.subject_recording import SubjectRecording

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
_SKLEARN_OK, _SKLEARN_ERR = True, ""
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold

default = Configuration()

def _build_pipeline(c_grid, inner_folds, scoring, seed):
    """Pipeline scikit-learn com TODA transformacao dentro do fold (evita vazamento)."""
    pipe = Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("sc", StandardScaler()),
            ("clf", LogisticRegression(penalty="l2", solver="liblinear", max_iter=5000))
        ])
    return GridSearchCV( pipe, 
                        {"clf__C": list(c_grid)},
                        cv=StratifiedKFold(inner_folds, shuffle=True, random_state=seed),
                        scoring=scoring, refit=True, n_jobs=1)

def aggregate_features(arrays: Dict[str, np.ndarray], channels: Sequence[str],
                       fcfg: Configuration) -> Dict[str, float]:
    """Agrega epoca->participante e canal->participante. p__ = confirmatorio, e__ = exploratorio."""
    agg = np.nanmedian if fcfg.aggregate == "median" else np.nanmean
    row: Dict[str, float] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in PRIMARY_FEATURES:
            v = arrays[k]
            row["p__" + k] = float(agg(np.nanmedian(v, axis=1)))
            for j, ch in enumerate(channels):
                row["e__%s_%s" % (k, ch)] = float(agg(v[:, j]))
            q = np.nanpercentile(np.nanmedian(v, axis=1), [25, 75])
            row["e__iqr_" + k] = float(q[1] - q[0])
        row["e__faa"] = float(agg(arrays["faa"]))
        row["e__aperiodic_r2"] = float(agg(np.nanmedian(arrays["aperiodic_r2"], axis=1)))
        row["e__alpha_detect_rate"] = float(np.mean(np.isfinite(arrays["alpha_peak_hz"])))
        for k in ("perm_entropy", "higuchi_fd", "c0_complexity"):
            v = arrays[k]
            row["e__" + k] = float(agg(np.nanmedian(v, axis=1)))
            for j, ch in enumerate(channels):
                row["e__%s_%s" % (k, ch)] = float(agg(v[:, j]))
    return row

def alpha_peak_from_flat(flat: np.ndarray, f_fit: np.ndarray,
                         fcfg: Configuration) -> np.ndarray:
    """Pico alfa detectado no espectro ACHATADO; sem pico identificavel -> NaN."""
    m = mask(f_fit, *fcfg.alpha_peak_band)
    if not m.any():
        return np.full(flat.shape[:-1], np.nan)
    sub = flat[..., m]
    idx = np.argmax(sub, axis=-1)
    height = np.take_along_axis(sub, idx[..., None], axis=-1)[..., 0]
    peak = f_fit[m][idx]
    return np.where(height >= fcfg.alpha_min_peak_log10, peak, np.nan)

def assert_no_circular_features(cols: Sequence[str]) -> None:
    """Escalas clinicas e o proprio rotulo nunca podem ser preditores."""
    bad = [c for c in cols if any(t in c.lower() for t in CLINICAL_SCALE_TOKENS)
           or c.lower() in {"label", "y", "type", "label_raw"}]
    if bad:
        raise RuntimeError("Features circulares/clinicas proibidas: %s" % bad[:8])

def assert_primary_only(cols: Sequence[str]) -> None:
    """Barreira dura confirmatorio/exploratorio."""
    bad = [c for c in cols if not c.startswith("p__")]
    if bad:
        raise RuntimeError("Colunas nao primarias no modelo confirmatorio: %s" % bad[:8])

def audit_expected_cohort(observed_ids: Sequence[str],
                          ccfg: Configuration) -> Dict[str, Any]:
    """B11: confere a coorte observada contra a lista canonica de IDs esperados.

    ``expected_n`` sozinho responde QUANTOS faltam; a lista de IDs responde QUAIS.
    Sem identidade declarada, um arquivo ausente e um arquivo nunca esperado sao
    indistinguiveis, e as contagens tipo CONSORT ficam sem denominador.
    """
    obs = sorted(set(map(str, observed_ids)))
    if not ccfg.expected_subject_ids:
        return {"expected_declared": False, "n_observed": len(obs)}
    exp = sorted(set(ccfg.expected_subject_ids))
    return {"expected_declared": True, "n_expected": len(exp), "n_observed": len(obs),
            "missing_ids": sorted(set(exp) - set(obs)),
            "unexpected_ids": sorted(set(obs) - set(exp)),
            "coverage": float(len(set(exp) & set(obs)) / len(exp))}

def band_share(
        frequencies: np.ndarray, 
        power_spectra: np.ndarray, 
        low_frequency: float, 
        high_frequency: float,
        reference: Tuple[float, float]
    ) -> np.ndarray:
    """Fracao adimensional da potencia de ``ref`` contida em [lo, hi)."""
    num = power_spectra[..., (frequencies >= low_frequency) & (frequencies < high_frequency)].sum(axis=-1)
    den = power_spectra[..., (frequencies >= reference[0]) & (frequencies < reference[1])].sum(axis=-1) + 1e-30
    return num / den

def band_share_median(
        frequencies:np.ndarray,
        power_spectra:np.ndarray,
        low_frequency:float,
        high_frequency: float,
        reference: Tuple[float,float]
                 ):
    return float(
        np.median(
            band_share(frequencies, power_spectra, low_frequency, high_frequency, reference)
            )
        )

def batch_confound_report(df: pd.DataFrame, acfg: Configuration) -> Dict[str, Any]:
    """Mede a colinearidade entre lote (prefixo de ID) e diagnostico, e o poder do lote.

    Se o prefixo separa perfeitamente os grupos, qualquer diferenca tecnica entre
    lotes e matematicamente indistinguivel do efeito clinico. O relatorio expoe
    esse fato e, quando V >= limiar, marca todos os resultados como confundidos.
    """
    if "id_prefix" not in df.columns:
        return {"available": False}
    v = cramers_v(df["id_prefix"], df["label"])
    tab = pd.crosstab(df["id_prefix"], df["label"]).to_dict()
    perfect = bool(df.groupby("id_prefix")["label"].nunique().max() == 1)
    return {"available": True, "cramers_v": v, "crosstab": tab,
            "prefix_perfectly_predicts_label": perfect,
            "batch_confounded": bool(perfect or v >= acfg.batch_collinearity_warn),
            "interpretation": ("Prefixo de ID determina o diagnostico: efeito de lote e "
                               "efeito clinico sao INSEPARAVEIS neste conjunto."
                               if perfect else
                               "Colinearidade parcial entre lote e diagnostico; ajustar e "
                               "reportar sensibilidade.")}

def bh_fdr(p: Sequence[float]) -> np.ndarray:
    """Correcao de Benjamini-Hochberg para multiplas comparacoes."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out

def bootstrap_ci(func, *arrays, n_boot: int = 2000, alpha: float = 0.05,
                 seed: int = 0) -> Tuple[float, float]:
    """IC percentil por reamostragem de PARTICIPANTES (unidade independente)."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        res = [a[rng.integers(0, len(a), len(a))] for a in arrays]
        v = func(*res)
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 100 * alpha / 2)),
            float(np.percentile(vals, 100 * (1 - alpha / 2))))

def c0_complexity(x: np.ndarray) -> float:
    """Complexidade C0: fracao de energia nao explicada pelas componentes dominantes."""
    X = np.fft.fft(x)
    p = np.abs(X) ** 2
    mean_p = p.mean()
    Y = np.where(p > mean_p, X, 0.0)
    y = np.fft.ifft(Y)
    num = float(np.sum(np.abs(x - y) ** 2))
    den = float(np.sum(np.abs(x) ** 2)) + 1e-30
    return num / den

def compute_psd(epochs: np.ndarray, fs: float, fcfg: Configuration):
    """PSD de Welch por epoca e canal."""
    nper = int(round(fcfg.welch_nperseg_seconds * fs))
    return sg.welch(epochs, fs=fs, nperseg=nper,
                  noverlap=int(round(nper * fcfg.welch_overlap)),
                  detrend="constant", axis=-1)

def confirmatory_group_tests(df: pd.DataFrame, acfg: Configuration) -> pd.DataFrame:
    """Testes da familia confirmatoria (p__), com efeito, IC e FDR de Benjamini-Hochberg."""
    cols = primary_columns(df)
    assert_primary_only(cols)
    assert_no_circular_features(cols)
    y = df["label"].to_numpy(int)
    rows = []
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        a, b = v[y == 1], v[y == 0]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if a.size < 3 or b.size < 3:
            rows.append({"feature": c, "n_mdd": a.size, "n_hc": b.size, "p_raw": np.nan,
                         "hedges_g": np.nan, "g_ci_low": np.nan, "g_ci_high": np.nan})
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        g = hedges_g(a, b)
        lo, hi = bootstrap_ci(hedges_g, a, b, n_boot=acfg.n_boot_effect,
                              seed=acfg.random_seed)
        rows.append({"feature": c, "n_mdd": int(a.size), "n_hc": int(b.size),
                     "median_mdd": float(np.median(a)), "median_hc": float(np.median(b)),
                     "u_statistic": float(u.statistic), "p_raw": float(u.pvalue),
                     "hedges_g": g, "g_ci_low": lo, "g_ci_high": hi})
    out = pd.DataFrame(rows)
    ok = out["p_raw"].notna()
    out.loc[ok, "p_fdr"] = bh_fdr(out.loc[ok, "p_raw"].to_numpy(float))
    out["significant_fdr"] = out["p_fdr"] < acfg.fdr_alpha
    return out.sort_values("p_raw").reset_index(drop=True)

def confound_adjusted_models(df: pd.DataFrame, acfg: Configuration) -> pd.DataFrame:
    """Regressao logistica por feature ajustada por idade, sexo e escolaridade."""
    import statsmodels.api as sm  # dependencia opcional
    cols = primary_columns(df)
    covs = [c for c in acfg.covariates if c in df.columns]
    rows = []
    for c in cols:
        sub = df[[c, "label"] + covs].apply(pd.to_numeric, errors="coerce").dropna()
        if len(sub) < 20:
            continue
        X = sm.add_constant(np.column_stack(
            [stats.zscore(sub[c].to_numpy(float))] +
            [stats.zscore(sub[k].to_numpy(float)) for k in covs]))
        try:
            res = sm.Logit(sub["label"].to_numpy(int), X).fit(disp=0)
            rows.append({"feature": c, "n": len(sub), "beta_z": float(res.params[1]),
                         "or_per_sd": float(np.exp(res.params[1])),
                         "p_value": float(res.pvalues[1]),
                         "ci_low": float(np.exp(res.conf_int()[1][0])),
                         "ci_high": float(np.exp(res.conf_int()[1][1])),
                         "covariates": ",".join(covs)})
        except Exception as exc:
            rows.append({"feature": c, "n": len(sub), "error": str(exc)[:120]})
    out = pd.DataFrame(rows)
    if "p_value" in out.columns and out["p_value"].notna().any():
        ok = out["p_value"].notna()
        out.loc[ok, "p_fdr"] = bh_fdr(out.loc[ok, "p_value"].to_numpy(float))
    return out

def canonical_subject_id(values) -> pd.Series:
    """Normaliza identificadores para 8 digitos com zeros a esquerda."""
    s = pd.Series(values).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(8)

def convert_numpy_to_pandas(data:np.ndarray, columns:Tuple[str]|None):
    return pd.DataFrame(np.transpose(data), columns=columns)

def cramers_v(x: Sequence, y: Sequence) -> float:
    """V de Cramer entre duas variaveis categoricas (0 = independentes, 1 = colineares)."""
    tab = pd.crosstab(pd.Series(x), pd.Series(y)).to_numpy()
    if tab.size == 0 or tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    n = tab.sum()
    return float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))

def discover_eeg_txt_files(eeg_dir: str | Path ) -> Tuple[List[Any], List[Dict[str, str]]]:
    """Descobre TXT validos; nomes fora do padrão vao para quarentena (nao abortam)."""
    root = Path(eeg_dir)
    if not root.is_dir():
        raise FileNotFoundError("Diretório com dados EEG inexistente: %s" % root)
    valid, quarantined = [], []
    for p in sorted(root.glob("*.txt")):
        try:
            parse_modma_filename(p.name)
            valid.append(p)
        except ValueError as exc:
            quarantined.append({"file": p.name, "stage": "filename", "error": str(exc)[:200]})
    if not valid:
        raise FileNotFoundError("Nenhum TXT MODMA valido encontrado.")
    return valid, quarantined

def epoch_feature_arrays(epochs: np.ndarray, fs: float, fcfg: Configuration,
                         channels: Sequence[str]) -> Dict[str, np.ndarray]:
    """Calcula todas as features no nivel EPOCA x CANAL."""
    f, psd = compute_psd(epochs, fs, fcfg)
    out: Dict[str, np.ndarray] = {}
    for name, val in relative_band_powers(f, psd, fcfg).items():
        out["rel_" + name] = val
    out["spec_entropy"] = spectral_entropy(f, psd, fcfg)
    slope, r2, flat, f_fit = robust_aperiodic_fit(f, psd, fcfg)
    out["aperiodic_slope"] = slope
    out["aperiodic_r2"] = r2
    out["alpha_peak_hz"] = alpha_peak_from_flat(flat, f_fit, fcfg)
    out["hjorth_mobility"], out["hjorth_complexity"] = hjorth_params(epochs)
    ne, nc, _ = epochs.shape
    pe = np.empty((ne, nc)); hf = np.empty((ne, nc)); c0 = np.empty((ne, nc))
    for i in range(ne):
        for j in range(nc):
            seg = epochs[i, j]
            pe[i, j] = permutation_entropy(seg, fcfg.perm_entropy_order, fcfg.perm_entropy_delay)
            hf[i, j] = higuchi_fd(seg, fcfg.higuchi_kmax)
            c0[i, j] = c0_complexity(seg)
    out["perm_entropy"], out["higuchi_fd"], out["c0_complexity"] = pe, hf, c0
    a = psd[..., mask(f, 8.0, 13.0)].sum(axis=-1)
    out["faa"] = np.log(a[:, 2] + 1e-30) - np.log(a[:, 0] + 1e-30)
    return out

def epoch_rejection_mask(epochs: np.ndarray, f: np.ndarray, psd: np.ndarray,
                         qc: Configuration, fcfg: Configuration) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Rejeicao de epoca relativa ao proprio sujeito, adimensional e cega ao rotulo."""
    finite = np.isfinite(epochs).all(axis=(1, 2))
    p2p = np.ptp(epochs, axis=2)
    jump = np.abs(np.diff(epochs, axis=2)).max(axis=2)
    nonflat = p2p.min(axis=1) > 0
    r_p2p = p2p / (np.median(p2p, axis=0) + 1e-30)
    r_jump = jump / (np.median(jump, axis=0) + 1e-30)
    lf = band_share(f, psd, fcfg.ocular_band[0], fcfg.ocular_band[1], fcfg.total_band)
    r_lf = lf / (np.median(lf, axis=0) + 1e-30)
    ok_hi = (r_p2p <= qc.p2p_ratio_max).all(axis=1)
    ok_lo = (r_p2p >= qc.p2p_ratio_min).all(axis=1)
    ok_jump = (r_jump <= qc.jump_ratio_max).all(axis=1)
    ok_lf = (r_lf <= qc.lf_share_ratio_max).all(axis=1)
    good = finite & nonflat & ok_hi & ok_lo & ok_jump & ok_lf
    diag = {"n_epochs": int(epochs.shape[0]), "n_good": int(good.sum()),
            "n_rejected": int((~good).sum()),
            "frac_rejected": float(1.0 - good.mean()) if good.size else float("nan"),
            "n_rej_nonfinite": int((~finite).sum()), "n_rej_flat": int((~nonflat).sum()),
            "n_rej_p2p_high": int((~ok_hi).sum()), "n_rej_p2p_low": int((~ok_lo).sum()),
            "n_rej_jump": int((~ok_jump).sum()), "n_rej_ocular_lf": int((~ok_lf).sum()),
            "rule": "relative_intrasubject_p2p_jump_lowfreq"}
    return good, diag

def epoch_signal(x: np.ndarray, fs: float, epoch_seconds: float,
                 overlap: float = 0.0) -> np.ndarray:
    """Segmenta a janela em epocas -> (n_epochs, n_channels, n_times)."""
    n = int(round(epoch_seconds * fs))
    step = max(1, int(round(n * (1.0 - overlap))))
    starts = list(range(0, x.shape[1] - n + 1, step))
    if not starts:
        raise ValueError("Nenhuma epoca gerada.")
    return np.stack([x[:, s:s + n] for s in starts], axis=0)

def extract_features_for_window(filtered: np.ndarray, rec: SubjectRecording, cfg: Configuration,
                                start_s: float, window_s: float) -> Dict[str, Any]:
    """Recorta a janela, epoca, aplica QC de epoca e agrega as features do participante."""
    # acq, pre, qc, fcfg = cfg.acquisition, cfg.preproc, cfg.qc, cfg.features
    win = extract_window(filtered, cfg.sampling_frequency, start_s, window_s)
    drift = within_window_drift(win, cfg.sampling_frequency, cfg.block_time)
    epochs = epoch_signal(win, cfg.sampling_frequency, cfg.epoch_seconds, cfg.epoch_overlap)
    f, psd = compute_psd(epochs, cfg.sampling_frequency, cfg)
    good, diag = epoch_rejection_mask(epochs, f, psd, cfg, cfg)
    if diag["frac_rejected"] > cfg.max_rejected_fraction:
        raise ValueError("%s: fracao de epocas rejeitadas %.3f > %.2f."
                         % (rec.subject_id, diag["frac_rejected"], cfg.max_rejected_fraction))
    if int(good.sum()) < cfg.min_good_epochs:
        raise ValueError("%s: %d epocas validas < %d."
                         % (rec.subject_id, int(good.sum()), cfg.min_good_epochs))
    arrays = epoch_feature_arrays(epochs[good], cfg.sampling_frequency, cfg, cfg.electrode_labels)
    row: Dict[str, Any] = {"subject_id": rec.subject_id, "source_name": rec.source_name,
                           "task": rec.task, "window_start_s": float(start_s),
                           "window_seconds": float(window_s),
                           "n_epochs_window": int(len(epochs)),
                           "n_epochs_used": int(good.sum()),
                           "frac_epochs_rejected": float(diag["frac_rejected"])}
    row.update({("qc_" + k): v for k, v in diag.items() if k.startswith("n_rej")})
    row.update({("drift_" + k.replace("drift_", "")): v for k, v in drift.items()})
    row.update(aggregate_features(arrays, cfg.electrode_labels, cfg))
    return row

def extract_window(filtered: np.ndarray, fs: float, start_s: float,
                   window_s: float) -> np.ndarray:
    """Recorta a janela [start_s, start_s + window_s) do sinal ja filtrado."""
    a, b = int(round(start_s * fs)), int(round((start_s + window_s) * fs))
    if filtered.shape[1] < b:
        raise ValueError("Sinal insuficiente: %d < %d amostras." % (filtered.shape[1], b))
    return np.ascontiguousarray(filtered[:, a:b])

def fingerprint(this) -> str:
    return hashlib.sha256(
        json.dumps(asdict(this), sort_keys=True, default=str).encode()).hexdigest()[:16]

def get_first_metadata(eeg_dir:str)->pd.DataFrame:
    """Carrega a primeira planilha xlsx encontrada no diretório, descartando colunas 'Unnamed'.
    """
    root = Path(eeg_dir)
    if not root.is_dir():
        raise FileNotFoundError("Diretório com dados EEG inexistente: %s" % root)
    files = sorted(root.glob("*.xlsx"))
    if(len(files)==0):
        raise FileNotFoundError("Nenhum arquivo xlsx encontrado: %s" % root)
    load = pd.read_excel(files[0])
    unnamed_columns = [col for col in load.columns if col.startswith('Unnamed')]
    for col in unnamed_columns:
        load.drop(col, axis=1, inplace=True) 
    return load, 

def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Tamanho de efeito padronizado com correcao para amostras pequenas."""
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    b = np.asarray(b, float); b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    na, nb = a.size, b.size
    sp = math.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp <= 0:
        return float("nan")
    d = (a.mean() - b.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (na + nb) - 9.0)
    return float(d * j)

def higuchi_fd(x: np.ndarray, kmax: int = 10) -> float:
    """Dimensao fractal de Higuchi: rugosidade da serie no tempo."""
    n = len(x)
    lk = []
    ks = []
    for k in range(1, kmax + 1):
        lm = []
        for m in range(k):
            idx = np.arange(m, n, k)
            if idx.size < 2:
                continue
            lm.append(np.abs(np.diff(x[idx])).sum() * (n - 1) / (idx.size - 1) / k)
        if lm:
            lk.append(np.log(np.mean(lm) + 1e-30)); ks.append(np.log(1.0 / k))
    if len(ks) < 2:
        return float("nan")
    return float(np.polyfit(ks, lk, 1)[0])

def hjorth_params(epochs: np.ndarray):
    """Mobilidade e complexidade de Hjorth: razoes de desvios, invariantes a escala."""
    d1 = np.diff(epochs, axis=-1)
    d2 = np.diff(d1, axis=-1)
    s0 = epochs.std(axis=-1) + 1e-30
    s1 = d1.std(axis=-1) + 1e-30
    s2 = d2.std(axis=-1) + 1e-30
    mob = s1 / s0
    return mob, (s2 / s1) / mob

def load_modma_metadata(path: str | Path) -> pd.DataFrame:
    """Le a planilha, canoniza IDs, mapeia MDD/HC e extrai covariaveis demograficas."""
    p = Path(path)
    raw = pd.read_csv(p) if p.suffix.lower() == ".csv" else pd.read_excel(p)
    if raw.empty:
        raise ValueError("Planilha de metadados vazia.")
    raw = raw.loc[:, [c for c in raw.columns if not str(c).lower().startswith("unnamed")]]
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]
    sidcols = [c for c in raw.columns
               if c in {"subject_id", "subject", "id", "participant_id"} or c.startswith("subject")]
    groupcols = [c for c in raw.columns
                 if c in {"type", "group", "class", "diagnosis", "diagnostic_group"}]
    if len(sidcols) != 1 or len(groupcols) != 1:
        raise ValueError("Esperada exatamente uma coluna de ID e uma de grupo.")
    out = raw.rename(columns={sidcols[0]: "subject_id", groupcols[0]: "label_raw"}).copy()
    out["subject_id"] = canonical_subject_id(out["subject_id"])
    if out["subject_id"].duplicated().any():
        raise ValueError("IDs duplicados nos metadados.")
    lab = out["label_raw"].astype(str).str.upper().str.strip()
    out["label"] = np.where(lab.str.startswith("MDD"), 1, np.where(lab.str.startswith("HC"), 0, -1))
    if (out["label"] < 0).any():
        raise ValueError("Rotulo nao mapeavel para MDD/HC.")
    for c in list(out.columns):
        if c.startswith("education"):
            out = out.rename(columns={c: "education_years"})
        elif c.startswith("gender") or c.startswith("sex"):
            out["sex"], sex_report = parse_sex_column(out[c])
    if "sex" in out.columns and not out["sex"].notna().any():
        warnings.warn("B9: coluna de sexo presente mas totalmente nao mapeada; os modelos "
                      "ajustados por confundidores ficariam silenciosamente vazios.",
                      RuntimeWarning)
    out["id_prefix"] = out["subject_id"].str[:4]
    keep = [c for c in ("subject_id", "label", "label_raw", "age", "sex",
                        "education_years", "id_prefix") if c in out.columns]
    return out[keep].drop_duplicates("subject_id").reset_index(drop=True)

def load_modma_txt(path: str | Path, cfg: Configuration, source_name: str) -> SubjectRecording:
    """Le um TXT MODMA de 3 canais e valida o esquema antes de devolver o registro."""
    
    sid, task = parse_modma_filename(source_name)
    arr = np.loadtxt(str(path), dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != cfg.expected_n_channels:
        raise ValueError("%s: esperado (N,%d), obtido %s."
                         % (source_name, cfg.expected_n_channels, arr.shape))
    if cfg.require_finite and not np.isfinite(arr).all():
        raise ValueError("%s: valores nao finitos." % source_name)
    if cfg.require_integer and not np.allclose(arr, np.round(arr)):
        raise ValueError("%s: valores nao inteiros." % source_name)
    if cfg.enforce_allowed_tasks and task not in cfg.allowed_tasks:
        raise ValueError("%s: tarefa '%s' fora do protocolo." % (source_name, task))
    arr, n_wrap = ft.fix_integer_wraparound(arr, cfg.container_bits)
    x = np.ascontiguousarray(arr.T)
    dur = x.shape[1] / cfg.sampling_frequency
    if dur > cfg.max_data_duration:
        raise ValueError("%s: duracao %.1f s acima do maximo." % (source_name, dur))
    rep = {"n_samples": int(x.shape[1]), "duration_s": float(dur),
           "n_wraparound_fixed": int(n_wrap),
           "wraparound_fraction": float(n_wrap / x.size),
           "dc_offset_counts": [float(v) for v in x.mean(axis=1)],
           "task_parsed": task}
    return SubjectRecording(sid, task, x, cfg.sampling_frequency, source_name, rep)

def mask(f: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return (f >= lo) & (f < hi)

def nested_cv_scores(X: np.ndarray, y: np.ndarray, acfg: Configuration,
                     repeats: Optional[int] = None) -> Dict[str, Any]:
    """Validacao cruzada ANINHADA no nivel do participante.

    O fold externo estima desempenho; o fold interno escolhe o hiperparametro.
    Nada e ajustado fora do fold, o que evita otimismo por vazamento de dados.
    """
    if not _SKLEARN_OK:
        raise RuntimeError("scikit-learn indisponivel: %s" % _SKLEARN_ERR)
    reps = acfg.outer_repeats if repeats is None else repeats
    cv = RepeatedStratifiedKFold(n_splits=acfg.outer_folds, n_repeats=reps,
                                 random_state=acfg.random_seed)
    oof = np.full((reps, len(y)), np.nan)
    for i, (tr, te) in enumerate(cv.split(X, y)):
        gs = _build_pipeline(acfg.logistic_c_grid, acfg.inner_folds,
                             acfg.inner_scoring, acfg.random_seed)
        gs.fit(X[tr], y[tr])
        oof[i // acfg.outer_folds, te] = gs.predict_proba(X[te])[:, 1]
    return {"oof_prob_per_repeat": oof, "prob_mean": np.nanmean(oof, axis=0)}

def nonstationarity_cv(
        data: np.ndarray, 
        sampling_frequency: float, 
        time_block: float = 10.0
    ) -> float:
    """Coeficiente de variacao da potencia em blocos: instabilidade do registro."""
    sample_size = int(round(time_block * sampling_frequency))
    if data.shape[1] < 3 * sample_size:
        return float("nan")
    nb = data.shape[1] // sample_size
    v = np.array([data[:, i * sample_size:(i + 1) * sample_size].var(axis=1).mean() for i in range(nb)])
    m = float(np.mean(v))
    return float(np.std(v) / m) if m > 0 else float("nan")

def parse_modma_filename(name: str) -> Tuple[str, str]:
    """Extrai (subject_id, task) do nome do arquivo; task ausente -> 'unspecified'."""
    m = NAME_REGEX.match(name.strip())
    if not m:
        raise ValueError("Nome fora do padrão MODMA: %r" % name)
    return m.group("sid"), (m.group("task") or 'undefined')

def parse_sex_column(values: pd.Series) -> Tuple[pd.Series, Dict[str, Any]]:
    """B9: parser de sexo robusto a codificacao numerica, com relatorio explicito.

    ERRO CORRIGIDO. A v23 fazia ``str[0]`` e mapeava apenas 'M'/'F'. Se a planilha
    codificasse sexo como 1/2 - convencao comum, inclusive no MODMA - o resultado
    era NaN para TODOS os participantes. Em ``confound_adjusted_models`` o
    ``dropna()`` esvaziava o dataframe, ``len(sub) < 20`` disparava e cada feature
    era pulada SEM erro: os modelos ajustados por confundidores simplesmente nao
    rodavam, e nada no relatorio deixava isso evidente.

    Convencao de saida: 1.0 = masculino, 0.0 = feminino, NaN = nao mapeado.
    A codificacao numerica assumida (1=M, 2=F) e declarada no relatorio para que
    possa ser contestada; se estiver invertida, o sinal do coeficiente de sexo
    inverte, sem afetar as demais covariaveis.
    """
    raw = pd.Series(values)
    s = raw.astype(str).str.strip().str.upper()
    out = pd.Series(np.nan, index=raw.index, dtype=float)
    out[s.str.startswith("M") | s.isin({"MALE", "MASCULINO", "H"})] = 1.0
    out[s.str.startswith("F") | s.isin({"FEMALE", "FEMININO", "W"})] = 0.0
    num = pd.to_numeric(raw, errors="coerce")
    coding = "text"
    unresolved = out.isna() & num.notna()
    if unresolved.any():
        vals = set(np.unique(num[unresolved].to_numpy()))
        if vals <= {1.0, 2.0}:
            out[unresolved] = np.where(num[unresolved] == 1.0, 1.0, 0.0); coding = "numeric_1M_2F"
        elif vals <= {0.0, 1.0}:
            out[unresolved] = num[unresolved].astype(float); coding = "numeric_1M_0F"
    report = {"n_total": int(len(raw)), "n_mapped": int(out.notna().sum()),
              "n_unmapped": int(out.isna().sum()), "coding_detected": coding,
              "distinct_input_values": sorted(map(str, pd.unique(raw.astype(str))))[:12]}
    if report["n_mapped"] == 0:
        warnings.warn("B9: nenhuma linha de sexo pode ser mapeada (valores: %s)."
                      % report["distinct_input_values"], RuntimeWarning)
    return out, report

def primary_columns(df: pd.DataFrame) -> List[str]:
    return sorted(c for c in df.columns if c.startswith("p__"))

def permutation_entropy(x: np.ndarray, order: int = 3, delay: int = 1) -> float:
    """Entropia de permutacao normalizada: previsibilidade dos padroes ordinais."""
    n = len(x) - (order - 1) * delay
    if n <= 1:
        return float("nan")
    idx = np.arange(order) * delay
    emb = np.stack([x[i + idx] for i in range(n)])
    perm = np.argsort(emb, axis=1)
    _, counts = np.unique(perm, axis=0, return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum() / np.log(math.factorial(order)))

def quality_control_metrics(
        raw_data:np.ndarray,
        filtered_data: np.ndarray, 
        cfg:Configuration = default,
        filtering_info={},
        subject_id = 'undefined',
        file_name = 'undefined',
        eeg_task = 'undefined',
        n_wraparound = 0,

        # fcfg: FeatureConfig, 
        # filt_diag: Dict[str, Any]
    ) -> Dict[str, Any]:
    """QC do REGISTRO INTEIRO: tecnica, adimensional e cega ao rotulo."""
    
    frequency_raw, power_raw = welch_power_spectral_density(raw_data - raw_data.mean(axis=1, keepdims=True), cfg.sampling_frequency)
    frequency_flt, power_flt = welch_power_spectral_density(filtered_data, cfg.sampling_frequency)

    low_24, high_24 = -(2 ** (cfg.bits_resolution - 1)), 2 ** (cfg.bits_resolution - 1) - 1
    diff = np.diff(raw_data, axis=1)

    correlation_matrix = np.corrcoef(filtered_data)
    offset_matrix = correlation_matrix[np.triu_indices(correlation_matrix.shape[0], k=1)]
    rms = float(np.median(np.sqrt((filtered_data ** 2).mean(axis=1))))
    return {
        "subject_id": subject_id, 
        "source_name": file_name, 
        "task": eeg_task,
        "n_samples": int(raw_data.shape[1]), 
        "duration_s": float(raw_data.shape[1] / cfg.sampling_frequency),
        "finite": bool(np.isfinite(raw_data).all()),
        "flat_channels": int(np.sum(np.ptp(raw_data, axis=1) == 0)),
        "n_wraparound_fixed": int(n_wraparound),
        "saturation_fraction": float((np.sum(raw_data <= low_24) + np.sum(raw_data >= high_24)) / raw_data.size),
        "zero_diff_fraction": float(np.mean(diff == 0)) if diff.size else float("nan"),
        "line_noise_ratio_raw": band_share_median(frequency_raw, power_raw, 49.0, 51.0, (1.0, 45.0)), # pq usou a referência de 1 a 45 e não de 1 a 40?
        "line_noise_ratio_post": float(filtering_info.get("line_ratio_post", np.nan)),
        "notch_applied": bool(filtering_info.get("notch_applied", False)),
        "ocular_index": band_share_median(frequency_flt, power_flt, *cfg.ocular_band, cfg.total_band),
        "muscle_ratio": band_share_median(frequency_flt, power_flt, *cfg.muscle_band, cfg.total_band),
        "nonstationarity_cv": nonstationarity_cv(filtered_data, cfg.sampling_frequency),
        "min_channel_corr": float(np.min(offset_matrix)) if offset_matrix.size else float("nan"),
        "rms_counts": rms, 
        "log_rms_counts": float(np.log10(rms + 1e-30))
    }

def relative_band_powers(f, psd, fcfg: Configuration) -> Dict[str, np.ndarray]:
    """Potencia RELATIVA por banda: adimensional, imune ao ganho de contato."""
    tot = psd[..., mask(f, *fcfg.total_band)].sum(axis=-1) + 1e-30
    return {n: psd[..., mask(f, lo, hi)].sum(axis=-1) / tot for n, lo, hi in fcfg.bands}

def rms_counts(data:np.ndarray):
    return float(np.median(np.sqrt((data ** 2).mean(axis=1))))

def robust_aperiodic_fit(f, psd, fcfg: Configuration):
    """Ajuste log-log robusto do fundo 1/f, removendo picos iterativamente.

    Devolve (inclinacao, R2, espectro achatado). A inclinacao e o expoente
    aperiodico com sinal invertido: quanto maior, mais "inclinado" o espectro.
    """
    m = mask(f, *fcfg.aperiodic_fit_band) & (f > 0)
    lf = np.log10(f[m])
    ly = np.log10(psd[..., m] + 1e-30)
    shape = ly.shape[:-1]
    lyf = ly.reshape(-1, ly.shape[-1])
    slopes = np.empty(lyf.shape[0]); r2s = np.empty(lyf.shape[0])
    flat = np.empty_like(lyf)
    X = np.vstack([lf, np.ones_like(lf)]).T
    for i in range(lyf.shape[0]):
        y = lyf[i]
        keep = np.ones_like(y, bool)
        beta = np.array([0.0, 0.0])
        for _ in range(fcfg.aperiodic_max_iter):
            beta, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            res = y - X @ beta
            sd = np.std(res[keep]) + 1e-30
            new = res <= fcfg.aperiodic_peak_sd * sd
            if new.sum() < 5 or np.array_equal(new, keep):
                keep = new if new.sum() >= 5 else keep
                break
            keep = new
        pred = X @ beta
        ss = np.sum((y[keep] - pred[keep]) ** 2)
        st = np.sum((y[keep] - y[keep].mean()) ** 2) + 1e-30
        slopes[i] = beta[0]; r2s[i] = 1.0 - ss / st
        flat[i] = y - pred
    return (slopes.reshape(shape), r2s.reshape(shape),
            flat.reshape(ly.shape), f[m])

def spectral_entropy(f, psd, fcfg: Configuration) -> np.ndarray:
    """Entropia de Shannon do espectro normalizado, em [0,1]: 1 = espectro plano."""
    sel = psd[..., mask(f, *fcfg.total_band)]
    p = sel / (sel.sum(axis=-1, keepdims=True) + 1e-30)
    return -(p * np.log(p + 1e-30)).sum(axis=-1) / np.log(p.shape[-1])

def sign_agreement_with_ci(g_ref: pd.Series, g_alt: pd.Series, n_boot: int = 2000,
                           seed: int = 0) -> Dict[str, Any]:
    """B7: concordancia de sinal entre janelas COM incerteza, sem veredito binario.

    ERRO CORRIGIDO. A v23 rotulava "consistente entre janelas" quando
    sign_agreement >= 0,7 sobre 9 features. Sob a hipotese nula de sinais
    aleatorios, P(>= 7 de 9 concordancias) = 0,090: o rotulo aparece por acaso em
    ~9% das vezes. Um limiar arbitrario sobre poucas features nao e verificacao de
    robustez. Agora reportamos a proporcao, o IC bootstrap por feature e o p-valor
    binomial exato contra o acaso (p = 0,5), sem rotulo qualitativo.
    """
    common = g_ref.index.intersection(g_alt.index)
    a = pd.to_numeric(g_ref[common], errors="coerce").to_numpy(float)
    b = pd.to_numeric(g_alt[common], errors="coerce").to_numpy(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = int(a.size)
    if n == 0:
        return {"n_features": 0, "sign_agreement": float("nan"),
                "ci95": [float("nan"), float("nan")], "p_binomial_vs_chance": float("nan")}
    agree = (np.sign(a) == np.sign(b)).astype(float)
    rng = np.random.default_rng(seed)
    boot = np.array([agree[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    k = int(agree.sum())
    p_bin = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)
    rho = stats.spearmanr(a, b).statistic if n >= 3 else float("nan")
    return {"n_features": n, "n_agreeing": k, "sign_agreement": float(agree.mean()),
            "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "p_binomial_vs_chance": p_bin,
            "spearman_effects": float(rho) if rho == rho else float("nan"),
            "note": ("Sem rotulo qualitativo: com poucas features a concordancia alta "
                     "e frequente por acaso (P(>=7/9) = 0,090).")}

def spectral_profile(
        frequencies: np.ndarray, 
        power_spectra:np.ndarray,
    ) -> np.ndarray:
    # """Vetor de participacao por banda (soma 1), medio entre canais: assinatura de estado."""
    # f, p = psd_continuous(seg, fs, nperseg_s=min(4.0, seg.shape[-1] / fs))

    v = np.array(
        [
            np.mean(
                band_share(frequencies, power_spectra, low, high, (1.0, 40.0))
            ) 
            for _, low, high in PROFILE_BANDS
        ])
    s = v.sum()
    return v / s if s > 0 else np.full(v.shape, 1.0 / v.size)

def spectral_profile_blocked(
        data: np.ndarray, 
        sampling_frequency: float|int, 
        time_window: float|int,
        cfg:Configuration = default
    ) -> pd.DataFrame:
    """
    Retorna o perfil espectral feito pelo método de welch 
    """
    block_size = int(round(time_window * sampling_frequency))
    steps = data.shape[1] // block_size
    rows = []
    for i in range(steps):
        segment = data[:, i * block_size:(i + 1) * block_size]
        frequencies, power_spectra = welch_power_spectral_density(segment, sampling_frequency, time_window=min(4.0, time_window))
        prof = spectral_profile(frequencies, power_spectra)
        rows.append({
            "block_index": i, 
            "t_start_s": i * time_window,
            "rms_counts": rms_counts(segment),
            "ocular_index": band_share_median(frequencies, power_spectra, *cfg.ocular_band,cfg.total_band),
            "muscle_ratio": band_share_median(frequencies, power_spectra, *cfg.muscle_band,cfg.total_band),
            **{"prof_%s" % nm: float(prof[j])
                        for j, (nm, _, _) in enumerate(PROFILE_BANDS)
                }
        })
    return pd.DataFrame(rows)

def verify_sampling_rate(x: np.ndarray, acq: Configuration) -> Dict[str, Any]:
    """B2: verificacao EMPIRICA da taxa de amostragem pelo pico de rede eletrica.

    ``fs`` nao e lida do arquivo: os TXT do MODMA sao matrizes de contas sem
    cabecalho. O descritor que a documenta especifica repouso de ~5 min, mas os
    registros tem ~1200 s - a mesma fonte e contraditada pelos dados. Portanto
    ``fs`` e SUPOSICAO, e esta funcao busca corroboracao independente.

    Metodo: sob a suposicao ``fs``, a interferencia de rede aparece em 50 Hz (China,
    MODMA) ou 60 Hz. Se o pico observado no espectro nao coincidir com nenhuma das
    frequencias candidatas dentro de ``tol``, a suposicao de ``fs`` fica SEM
    corroboracao e todo o eixo de frequencia escala por um fator desconhecido.

    LIMITACAO: a ausencia de pico de rede NAO refuta ``fs`` (o aparelho pode ter
    filtro de rede interno); apenas deixa a suposicao sem verificacao. A presenca
    do pico na posicao esperada e evidencia consistente, nao prova.
    """
    y = np.asarray(x, float)
    y = y - y.mean(axis=-1, keepdims=True)
    f, p = welch_power_spectral_density(y, acq.sampling_frequency, time_window=8.0)
    p = np.median(np.atleast_2d(p), axis=0)
    band = (f >= 35.0) & (f <= min(70.0, acq.sampling_frequency / 2.0 - 1.0))
    if band.sum() < 8:
        return {"fs_assumed_hz": float(acq.sampling_frequency), "fs_verification": "insufficient_band",
                "fs_supported": False}
    fb, pb = f[band], np.log10(p[band] + 1e-30)
    
    i = int(np.argmax(pb))
    peak_hz = float(fb[i])
    prominence = float(pb[i] - np.median(pb))
    cand = [c for c in acq.fs_verification_line_hz
            if abs(peak_hz - c) <= acq.fs_verification_tol_hz]
    has_peak = prominence >= acq.fs_verification_min_prominence_log10
    if has_peak and cand:
        status, supported = "line_peak_at_expected_frequency", True
    elif has_peak:
        status, supported = "line_peak_at_unexpected_frequency", False
    else:
        status, supported = "no_line_peak_detected_inconclusive", False
    implied = [float(acq.sampling_frequency * c / peak_hz) for c in acq.fs_verification_line_hz] if has_peak else []
    return {"implied_fs_if_peak_is_line_hz": implied,"fs_assumed_hz": float(acq.sampling_frequency), "fs_is_assumption": bool(acq.fs_is_assumption),
            "line_peak_hz_under_assumed_fs": peak_hz,
            "line_peak_prominence_log10": prominence,
            "fs_verification": status, "fs_supported": bool(supported)
            }

def welch_power_spectral_density(
        data: np.ndarray, 
        sampling_frequency: float|int, 
        time_window: float|int = 4.0):
    """PSD pelo método de Welch do sinal continuo (janelas de nperseg_s segundos)."""
    segment_length = min(int(round(time_window * sampling_frequency)), data.shape[-1])
    return sg.welch(data, fs=sampling_frequency, nperseg=segment_length, noverlap=segment_length // 2, detrend="constant", axis=-1)

def within_window_drift(window: np.ndarray, fs: float,
                        block_seconds: float) -> Dict[str, float]:
    """A5: mede se o ESTADO muda dentro da janela escolhida (deriva residual)."""
    prof = spectral_profile_blocked(window, fs, min(block_seconds, window.shape[1] / fs / 3.0))
    if len(prof) < 2:
        return {"drift_log2_rms": float("nan"), "drift_js": float("nan"),
                "drift_rms_spearman": float("nan")}
    M = an.profile_matrix(prof)
    med = np.median(M, axis=0)
    r = prof["rms_counts"].to_numpy(float)
    return {"drift_log2_rms": float(np.log2((r.max() + 1e-30) / (r.min() + 1e-30))),
            "drift_js": float(max(an.jensen_shannon(M[i], med) for i in range(len(prof)))),
            "drift_rms_spearman": float(stats.spearmanr(prof["t_start_s"], r).statistic)}

