from acquisition import AcquisitionConfig
from temporal_protocol import TemporalProtocolConfig
from pre_processor import PreprocConfig
from cohort import CohortConfig
from feature import FeatureConfig
from qc import QCConfig
from analysis import AnalysisConfig
from schema import SchemaConfig
from subject_recording import SubjectRecording
from run_config import RunConfig
from configuration import Configuration

from vars import TASK_UNSPECIFIED, NAME_REGEX, PROFILE_BANDS, PRIMARY_FEATURES, CLINICAL_SCALE_TOKENS

from typing import Tuple, Sequence, Dict, Any, Optional, List

from scipy import stats
from scipy.signal import butter, filtfilt, iirnotch, welch as _welch
import numpy as np
import pandas as pd
import warnings
import math
import sys
import hashlib
import platform
from datetime import datetime, timezone

from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import (average_precision_score, balanced_accuracy_score, brier_score_loss, confusion_matrix, f1_score, matthews_corrcoef, roc_auc_score)

_SKLEARN_OK, _SKLEARN_ERR = True, ""


def discover_eeg_files(eeg_dir: str | Path) -> Tuple[List[Path], List[Dict[str, str]]]:
    """Descobre TXT validos; nomes fora do padrao vao para quarentena (nao abortam)."""
    root = Path(eeg_dir)
    if not root.is_dir():
        raise FileNotFoundError("Diretorio com dados EEG inexistente: %s" % root)
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

def parse_modma_filename(name: str) -> Tuple[str, str]:
    """Extrai (subject_id, task) do nome do arquivo; task ausente -> 'unspecified'."""
    m = NAME_REGEX.match(name.strip())
    if not m:
        raise ValueError("Nome fora do padrao MODMA: %r" % name)
    return m.group("sid"), (m.group("task") or TASK_UNSPECIFIED)


def jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Divergencia de Jensen-Shannon (base 2) entre dois perfis espectrais em [0,1]."""
    p = np.clip(np.asarray(p, float), 1e-12, None); p = p / p.sum()
    q = np.clip(np.asarray(q, float), 1e-12, None); q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: float(np.sum(a * np.log2(a / b)))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def cramers_v(x: Sequence, y: Sequence) -> float:
    """V de Cramer entre duas variaveis categoricas (0 = independentes, 1 = colineares)."""
    tab = pd.crosstab(pd.Series(x), pd.Series(y)).to_numpy()
    if tab.size == 0 or tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    n = tab.sum()
    return float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))


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


def verify_sampling_rate(x: np.ndarray, acq: AcquisitionConfig|Configuration) -> Dict[str, Any]:
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
    f, p = psd_continuous(y, acq.fs, nperseg_s=8.0)
    p = np.median(np.atleast_2d(p), axis=0)
    band = (f >= 35.0) & (f <= min(70.0, acq.fs / 2.0 - 1.0))
    if band.sum() < 8:
        return {"fs_assumed_hz": float(acq.fs), "fs_verification": "insufficient_band",
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
    implied = [float(acq.fs * c / peak_hz) for c in acq.fs_verification_line_hz] if has_peak else []
    return {"fs_assumed_hz": float(acq.fs), "fs_is_assumption": bool(acq.fs_is_assumption),
            "line_peak_hz_under_assumed_fs": peak_hz,
            "line_peak_prominence_log10": prominence,
            "fs_verification": status, "fs_supported": bool(supported),
            "implied_fs_if_peak_is_line_hz": implied}


def psd_continuous(x: np.ndarray, fs: float, nperseg_s: float = 4.0):
    """PSD de Welch do sinal continuo (janelas de nperseg_s segundos)."""
    nper = min(int(round(nperseg_s * fs)), x.shape[-1])
    return _welch(x, fs=fs, nperseg=nper, noverlap=nper // 2, detrend="constant", axis=-1)



def preprocess_continuous(x_counts: np.ndarray, acq: AcquisitionConfig,
                          pre: PreprocConfig) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Notch condicional + passa-banda de fase zero + corte de borda, no sinal continuo."""
    y = np.asarray(x_counts, float)
    nyq = acq.fs / 2.0
    pre_ratio = line_noise_ratio(y - y.mean(axis=-1, keepdims=True), acq.fs, pre.notch_hz or 50.0)
    if pre.notch_mode == "always":
        do_notch = pre.notch_hz is not None
    elif pre.notch_mode == "never":
        do_notch = False
    else:
        do_notch = (pre.notch_hz is not None) and (pre_ratio > pre.notch_line_ratio_threshold)
    if do_notch:
        bn, an = iirnotch(pre.notch_hz, Q=pre.notch_q, fs=acq.fs)
        y = filtfilt(bn, an, y, axis=-1)
    lo, hi = pre.bandpass_hz
    b, a = butter(pre.filter_order, [lo / nyq, hi / nyq], btype="band")
    y = filtfilt(b, a, y, axis=-1)
    k = int(round(pre.edge_trim_seconds * acq.fs))
    if k > 0:
        if y.shape[1] <= 2 * k:
            raise ValueError("Sinal curto demais para o corte de borda.")
        y = y[:, k:-k]
    y = np.ascontiguousarray(y - y.mean(axis=1, keepdims=True))
    diag = {"notch_mode": pre.notch_mode, "notch_applied": bool(do_notch),
            "line_ratio_pre": float(pre_ratio),
            "line_ratio_post": float(line_noise_ratio(y, acq.fs, pre.notch_hz or 50.0)),
            "edge_trim_s": float(pre.edge_trim_seconds),
            "n_samples_after_trim": int(y.shape[1])}
    return y, diag


"""Gera sinal sintetico de 3 canais com transitorio inicial de amplitude."""
def synth_recording(fs: float, seconds: float, alpha_hz: float = 10.0,
                     settle_s: float = 120.0, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = int(fs * seconds)
    t = np.arange(n) / fs
    env = 1.0 + 3.0 * np.exp(-t / max(settle_s / 3.0, 1.0))  # acomodacao decrescente
    base = (np.sin(2 * np.pi * alpha_hz * t) + 0.5 * rng.standard_normal(n))
    x = np.stack([env * (base + 0.1 * rng.standard_normal(n)) for _ in range(3)])
    return np.round(x * 1000.0)


def block_profile(filtered: np.ndarray, fs: float, block_seconds: float) -> pd.DataFrame:
    """A5: descreve o registro em blocos consecutivos (RMS, ocular, bandas).

    Cada linha e um bloco; as colunas descrevem o ESTADO do sinal naquele bloco.
    E a materia-prima para decidir quando o registro se estabiliza.
    """
    n = int(round(block_seconds * fs))
    nb = filtered.shape[1] // n
    rows = []
    for i in range(nb):
        seg = filtered[:, i * n:(i + 1) * n]
        f, p = psd_continuous(seg, fs, nperseg_s=min(4.0, block_seconds))
        prof = spectral_profile(seg, fs)
        rows.append({"block_index": i, "t_start_s": i * block_seconds,
                     "rms_counts": float(np.median(np.sqrt((seg ** 2).mean(axis=1)))),
                     "ocular_index": float(np.median(band_share(f, p, 0.5, 3.0, (1.0, 40.0)))),
                     "muscle_ratio": float(np.median(band_share(f, p, 20.0, 40.0, (1.0, 40.0)))),
                     **{"prof_%s" % nm: float(prof[j])
                        for j, (nm, _, _) in enumerate(PROFILE_BANDS)}})
    return pd.DataFrame(rows)

def cohort_terminal_profile(block_frames: Sequence[pd.DataFrame],
                            tcfg: TemporalProtocolConfig) -> Optional[np.ndarray]:
    """B3: perfil espectral de referencia EXTERNO, mediana da coorte nos blocos finais.

    Calculado uma unica vez, sobre todos os sujeitos, sem qualquer acesso ao rotulo.
    E etapa de COORTE (nao pre-processamento individual) e deve ser descrita como tal.
    """
    rows = []
    for prof in block_frames:
        if len(prof) < 4:
            continue
        k = max(2, min(int(tcfg.n_terminal_blocks), len(prof) // 2))
        rows.append(np.median(profile_matrix(prof.iloc[-k:]), axis=0))
    if not rows:
        return None
    ref = np.median(np.vstack(rows), axis=0)
    total = ref.sum()
    return ref / total if total > 0 else ref

def subject_settling_time(prof: pd.DataFrame, tcfg: TemporalProtocolConfig, cohort_ref_profile: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """B3: instante a partir do qual o registro do sujeito fica ESTAVEL.

    PROBLEMA CORRIGIDO. Na v23 a referencia de estabilidade era a mediana da
    SEGUNDA METADE DO PROPRIO REGISTRO. Sob deriva monotona - e 31 dos 55
    registros tem Spearman(RMS, tempo) < -0,5 - a referencia tambem deriva: o
    metodo compara o sinal com um alvo movel e nao pode, por construcao, detectar
    deriva que atravessa o registro inteiro. Isso tende a declarar estabilidade
    cedo demais.

    CORRECAO. A referencia ESPECTRAL passa a ser EXTERNA ao sujeito: mediana da
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
                "ref_mode": tcfg.settling_ref_mode}
    k = max(2, min(int(tcfg.n_terminal_blocks), n // 2))
    tail = prof.iloc[-k:]
    rms_ref = float(np.median(tail["rms_counts"])) + 1e-30
    M = profile_matrix(prof)
    if tcfg.settling_ref_mode == "cohort" and cohort_ref_profile is not None:
        prof_ref = np.asarray(cohort_ref_profile, float)
        ref_used = "cohort_terminal_blocks"
    else:
        prof_ref = np.median(profile_matrix(tail), axis=0)
        ref_used = "self_terminal_blocks"
    dev_log2 = np.abs(np.log2((prof["rms_counts"].to_numpy(float) + 1e-30) / rms_ref))
    dev_js = np.array([jensen_shannon(M[i], prof_ref) for i in range(n)])
    stable = (dev_log2 <= tcfg.settle_tol_log2) & (dev_js <= tcfg.settle_tol_js)
    idx = len(stable)
    for i in range(len(stable) - 1, -1, -1):
        if stable[i]:
            idx = i
        else:
            break
    t = float(prof["t_start_s"].iloc[idx]) if idx < len(prof) else float("nan")
    tail_rho = float(stats.spearmanr(tail["t_start_s"], tail["rms_counts"]).statistic) \
        if k >= 3 else float("nan")
    spans = bool(np.isfinite(tail_rho) and abs(tail_rho) >= tcfg.tail_trend_warn)
    return {"settling_time_s": t, "n_blocks": int(n),
            "frac_blocks_stable": float(np.mean(stable)),
            "rms_trend_spearman": float(stats.spearmanr(prof["t_start_s"],
                                                        prof["rms_counts"]).statistic),
            "tail_rms_trend_spearman": tail_rho,
            "drift_spans_record": spans,
            "ref_mode": tcfg.settling_ref_mode, "reference_used": ref_used,
            "n_terminal_blocks": int(k),
            "dev_log2_first_block": float(dev_log2[0]), "dev_js_first_block": float(dev_js[0]),
            "settling_status": ("ok" if np.isfinite(t) else "never_stable")}

def line_noise_ratio(x: np.ndarray, fs: float, hz: float = 50.0, bw: float = 1.0,
                     ref: Tuple[float, float] = (1.0, 45.0)) -> float:
    """Evidencia MEDIDA de interferencia de rede eletrica; base da decisao de notch."""
    f, p = psd_continuous(x, fs)
    return float(np.median(np.atleast_1d(band_share(f, p, hz - bw, hz + bw, ref))))

def band_share(f: np.ndarray, psd: np.ndarray, lo: float, hi: float,
               ref: Tuple[float, float]) -> np.ndarray:
    """Fracao adimensional da potencia de ``ref`` contida em [lo, hi)."""
    num = psd[..., (f >= lo) & (f < hi)].sum(axis=-1)
    den = psd[..., (f >= ref[0]) & (f < ref[1])].sum(axis=-1) + 1e-30
    return num / den

def spectral_profile(seg: np.ndarray, fs: float) -> np.ndarray:
    """Vetor de participacao por banda (soma 1), medio entre canais: assinatura de estado."""
    f, p = psd_continuous(seg, fs, nperseg_s=min(4.0, seg.shape[-1] / fs))
    v = np.array([np.mean(band_share(f, p, lo, hi, (1.0, 40.0))) for _, lo, hi in PROFILE_BANDS])
    s = v.sum()
    return v / s if s > 0 else np.full(v.shape, 1.0 / v.size)

def profile_matrix(df: pd.DataFrame) -> np.ndarray:
    return df[["prof_%s" % nm for nm, _, _ in PROFILE_BANDS]].to_numpy(float)


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


def audit_expected_cohort(observed_ids: Sequence[str],
                          ccfg: CohortConfig) -> Dict[str, Any]:
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

def epoch_signal(x: np.ndarray, fs: float, epoch_seconds: float,
                 overlap: float = 0.0) -> np.ndarray:
    """Segmenta a janela em epocas -> (n_epochs, n_channels, n_times)."""
    n = int(round(epoch_seconds * fs))
    step = max(1, int(round(n * (1.0 - overlap))))
    starts = list(range(0, x.shape[1] - n + 1, step))
    if not starts:
        raise ValueError("Nenhuma epoca gerada.")
    return np.stack([x[:, s:s + n] for s in starts], axis=0)

def extract_window(filtered: np.ndarray, fs: float, start_s: float,
                   window_s: float) -> np.ndarray:
    """Recorta a janela [start_s, start_s + window_s) do sinal ja filtrado."""
    a, b = int(round(start_s * fs)), int(round((start_s + window_s) * fs))
    if filtered.shape[1] < b:
        raise ValueError("Sinal insuficiente: %d < %d amostras." % (filtered.shape[1], b))
    return np.ascontiguousarray(filtered[:, a:b])

def compute_psd(epochs: np.ndarray, fs: float, fcfg: FeatureConfig):
    """PSD de Welch por epoca e canal."""
    nper = int(round(fcfg.welch_nperseg_seconds * fs))
    return _welch(epochs, fs=fs, nperseg=nper,
                  noverlap=int(round(nper * fcfg.welch_overlap)),
                  detrend="constant", axis=-1)


def epoch_rejection_mask(epochs: np.ndarray, f: np.ndarray, psd: np.ndarray,
                         qc: QCConfig, fcfg: FeatureConfig) -> Tuple[np.ndarray, Dict[str, Any]]:
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

def epoch_feature_arrays(epochs: np.ndarray, fs: float, fcfg: FeatureConfig,
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

def aggregate_features(arrays: Dict[str, np.ndarray], channels: Sequence[str],
                       fcfg: FeatureConfig) -> Dict[str, float]:
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


def relative_band_powers(f, psd, fcfg: FeatureConfig) -> Dict[str, np.ndarray]:
    """Potencia RELATIVA por banda: adimensional, imune ao ganho de contato."""
    tot = psd[..., mask(f, *fcfg.total_band)].sum(axis=-1) + 1e-30
    return {n: psd[..., mask(f, lo, hi)].sum(axis=-1) / tot for n, lo, hi in fcfg.bands}


def mask(f: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return (f >= lo) & (f < hi)


def spectral_entropy(f, psd, fcfg: FeatureConfig) -> np.ndarray:
    """Entropia de Shannon do espectro normalizado, em [0,1]: 1 = espectro plano."""
    sel = psd[..., mask(f, *fcfg.total_band)]
    p = sel / (sel.sum(axis=-1, keepdims=True) + 1e-30)
    return -(p * np.log(p + 1e-30)).sum(axis=-1) / np.log(p.shape[-1])



def robust_aperiodic_fit(f, psd, fcfg: FeatureConfig):
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



def alpha_peak_from_flat(flat: np.ndarray, f_fit: np.ndarray,
                         fcfg: FeatureConfig) -> np.ndarray:
    """Pico alfa detectado no espectro ACHATADO; sem pico identificavel -> NaN."""
    m = mask(f_fit, *fcfg.alpha_peak_band)
    if not m.any():
        return np.full(flat.shape[:-1], np.nan)
    sub = flat[..., m]
    idx = np.argmax(sub, axis=-1)
    height = np.take_along_axis(sub, idx[..., None], axis=-1)[..., 0]
    peak = f_fit[m][idx]
    return np.where(height >= fcfg.alpha_min_peak_log10, peak, np.nan)


def hjorth_params(epochs: np.ndarray):
    """Mobilidade e complexidade de Hjorth: razoes de desvios, invariantes a escala."""
    d1 = np.diff(epochs, axis=-1)
    d2 = np.diff(d1, axis=-1)
    s0 = epochs.std(axis=-1) + 1e-30
    s1 = d1.std(axis=-1) + 1e-30
    s2 = d2.std(axis=-1) + 1e-30
    mob = s1 / s0
    return mob, (s2 / s1) / mob


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


def assert_primary_only(cols: Sequence[str]) -> None:
    """Barreira dura confirmatorio/exploratorio."""
    bad = [c for c in cols if not c.startswith("p__")]
    if bad:
        raise RuntimeError("Colunas nao primarias no modelo confirmatorio: %s" % bad[:8])


def assert_no_circular_features(cols: Sequence[str]) -> None:
    """Escalas clinicas e o proprio rotulo nunca podem ser preditores."""
    bad = [c for c in cols if any(t in c.lower() for t in CLINICAL_SCALE_TOKENS)
           or c.lower() in {"label", "y", "type", "label_raw"}]
    if bad:
        raise RuntimeError("Features circulares/clinicas proibidas: %s" % bad[:8])

def synth_recording(fs: float, seconds: float, alpha_hz: float = 10.0,
                     settle_s: float = 120.0, seed: int = 0) -> np.ndarray:
    """Gera sinal sintetico de 3 canais com transitorio inicial de amplitude."""
    rng = np.random.default_rng(seed)
    n = int(fs * seconds)
    t = np.arange(n) / fs
    env = 1.0 + 3.0 * np.exp(-t / max(settle_s / 3.0, 1.0))  # acomodacao decrescente
    base = (np.sin(2 * np.pi * alpha_hz * t) + 0.5 * rng.standard_normal(n))
    x = np.stack([env * (base + 0.1 * rng.standard_normal(n)) for _ in range(3)])
    return np.round(x * 1000.0)

def canonical_subject_id(values) -> pd.Series:
    """Normaliza identificadores para 8 digitos com zeros a esquerda."""
    s = pd.Series(values).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(8)


def load_modma_txt(path: str | Path, acq: AcquisitionConfig, source_name: str,
                   schema: SchemaConfig) -> SubjectRecording:
    """Le um TXT MODMA de 3 canais e valida o esquema antes de devolver o registro."""
    sid, task = parse_modma_filename(source_name)
    arr = np.loadtxt(str(path), dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != schema.expected_n_channels:
        raise ValueError("%s: esperado (N,%d), obtido %s."
                         % (source_name, schema.expected_n_channels, arr.shape))
    if schema.require_finite and not np.isfinite(arr).all():
        raise ValueError("%s: valores nao finitos." % source_name)
    if schema.require_integer and not np.allclose(arr, np.round(arr)):
        raise ValueError("%s: valores nao inteiros." % source_name)
    if schema.enforce_allowed_tasks and task not in schema.allowed_tasks:
        raise ValueError("%s: tarefa '%s' fora do protocolo." % (source_name, task))
    arr, n_wrap = fix_integer_wraparound(arr, acq.container_bits)
    x = np.ascontiguousarray(arr.T)
    dur = x.shape[1] / acq.fs
    if dur > schema.max_duration_s:
        raise ValueError("%s: duracao %.1f s acima do maximo." % (source_name, dur))
    rep = {"n_samples": int(x.shape[1]), "duration_s": float(dur),
           "n_wraparound_fixed": int(n_wrap),
           "wraparound_fraction": float(n_wrap / x.size),
           "dc_offset_counts": [float(v) for v in x.mean(axis=1)],
           "task_parsed": task}
    return SubjectRecording(sid, task, x, acq.fs, source_name, rep)


def subject_qc_metrics(rec: SubjectRecording, filtered: np.ndarray, acq: AcquisitionConfig,
                       fcfg: FeatureConfig, filt_diag: Dict[str, Any]) -> Dict[str, Any]:
    """QC do REGISTRO INTEIRO: tecnica, adimensional e cega ao rotulo."""
    raw = np.asarray(rec.data_counts, float)
    f_raw, p_raw = psd_continuous(raw - raw.mean(axis=1, keepdims=True), acq.fs)
    f_flt, p_flt = psd_continuous(filtered, acq.fs)
    lo24, hi24 = -(2 ** (acq.adc_bits - 1)), 2 ** (acq.adc_bits - 1) - 1
    diff = np.diff(raw, axis=1)
    corr = np.corrcoef(filtered)
    off = corr[np.triu_indices(corr.shape[0], k=1)]
    rms = float(np.median(np.sqrt((filtered ** 2).mean(axis=1))))
    return {"subject_id": rec.subject_id, "source_name": rec.source_name, "task": rec.task,
            "n_samples": int(raw.shape[1]), "duration_s": float(raw.shape[1] / acq.fs),
            "finite": bool(np.isfinite(raw).all()),
            "flat_channels": int(np.sum(np.ptp(raw, axis=1) == 0)),
            "n_wraparound_fixed": int(rec.parse_report.get("n_wraparound_fixed", 0)),
            "saturation_fraction": float((np.sum(raw <= lo24) + np.sum(raw >= hi24)) / raw.size),
            "zero_diff_fraction": float(np.mean(diff == 0)) if diff.size else float("nan"),
            "line_noise_ratio_raw": float(np.median(band_share(f_raw, p_raw, 49.0, 51.0, (1.0, 45.0)))),
            "line_noise_ratio_post": float(filt_diag.get("line_ratio_post", np.nan)),
            "notch_applied": bool(filt_diag.get("notch_applied", False)),
            "ocular_index": float(np.median(band_share(f_flt, p_flt, *fcfg.ocular_band,
                                                       ref=fcfg.total_band))),
            "muscle_ratio": float(np.median(band_share(f_flt, p_flt, *fcfg.muscle_band,
                                                       ref=fcfg.total_band))),
            "nonstationarity_cv": nonstationarity_cv(filtered, acq.fs),
            "min_channel_corr": float(np.min(off)) if off.size else float("nan"),
            "rms_counts": rms, "log_rms_counts": float(np.log10(rms + 1e-30))}

def nonstationarity_cv(x: np.ndarray, fs: float, block_s: float = 10.0) -> float:
    """Coeficiente de variacao da potencia em blocos: instabilidade do registro."""
    n = int(round(block_s * fs))
    if x.shape[1] < 3 * n:
        return float("nan")
    nb = x.shape[1] // n
    v = np.array([x[:, i * n:(i + 1) * n].var(axis=1).mean() for i in range(nb)])
    m = float(np.mean(v))
    return float(np.std(v) / m) if m > 0 else float("nan")



def extract_features_for_window(filtered: np.ndarray, rec: SubjectRecording, cfg: RunConfig,
                                start_s: float, window_s: float) -> Dict[str, Any]:
    """Recorta a janela, epoca, aplica QC de epoca e agrega as features do participante."""
    acq, pre, qc, fcfg = cfg.acquisition, cfg.preproc, cfg.qc, cfg.features
    win = extract_window(filtered, acq.fs, start_s, window_s)
    drift = within_window_drift(win, acq.fs, cfg.temporal.block_seconds)
    epochs = epoch_signal(win, acq.fs, pre.epoch_seconds, pre.epoch_overlap)
    f, psd = compute_psd(epochs, acq.fs, fcfg)
    good, diag = epoch_rejection_mask(epochs, f, psd, qc, fcfg)
    if diag["frac_rejected"] > qc.max_rejected_fraction:
        raise ValueError("%s: fracao de epocas rejeitadas %.3f > %.2f."
                         % (rec.subject_id, diag["frac_rejected"], qc.max_rejected_fraction))
    if int(good.sum()) < fcfg.min_good_epochs:
        raise ValueError("%s: %d epocas validas < %d."
                         % (rec.subject_id, int(good.sum()), fcfg.min_good_epochs))
    arrays = epoch_feature_arrays(epochs[good], acq.fs, fcfg, acq.channel_names)
    row: Dict[str, Any] = {"subject_id": rec.subject_id, "source_name": rec.source_name,
                           "task": rec.task, "window_start_s": float(start_s),
                           "window_seconds": float(window_s),
                           "n_epochs_window": int(len(epochs)),
                           "n_epochs_used": int(good.sum()),
                           "frac_epochs_rejected": float(diag["frac_rejected"])}
    row.update({("qc_" + k): v for k, v in diag.items() if k.startswith("n_rej")})
    row.update({("drift_" + k.replace("drift_", "")): v for k, v in drift.items()})
    row.update(aggregate_features(arrays, acq.channel_names, cfg.features))
    return row

def within_window_drift(window: np.ndarray, fs: float,
                        block_seconds: float) -> Dict[str, float]:
    """A5: mede se o ESTADO muda dentro da janela escolhida (deriva residual)."""
    prof = block_profile(window, fs, min(block_seconds, window.shape[1] / fs / 3.0))
    if len(prof) < 2:
        return {"drift_log2_rms": float("nan"), "drift_js": float("nan"),
                "drift_rms_spearman": float("nan")}
    M = profile_matrix(prof)
    med = np.median(M, axis=0)
    r = prof["rms_counts"].to_numpy(float)
    return {"drift_log2_rms": float(np.log2((r.max() + 1e-30) / (r.min() + 1e-30))),
            "drift_js": float(max(jensen_shannon(M[i], med) for i in range(len(prof)))),
            "drift_rms_spearman": float(stats.spearmanr(prof["t_start_s"], r).statistic)}

def audit_expected_cohort(observed_ids: Sequence[str],
                          ccfg: CohortConfig) -> Dict[str, Any]:
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

def confirmatory_group_tests(df: pd.DataFrame, acfg: AnalysisConfig) -> pd.DataFrame:
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


def confound_adjusted_models(df: pd.DataFrame, acfg: AnalysisConfig) -> pd.DataFrame:
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


def batch_confound_report(df: pd.DataFrame, acfg: AnalysisConfig) -> Dict[str, Any]:
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


def primary_columns(df: pd.DataFrame) -> List[str]:
    return sorted(c for c in df.columns if c.startswith("p__"))




def nested_cv_scores(X: np.ndarray, y: np.ndarray, acfg: AnalysisConfig,
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


def metrics_from_probs(y: np.ndarray, p: np.ndarray, acfg: AnalysisConfig,
                       oof_per_repeat: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Metricas de classificacao com IC bootstrap por participante.

    B11: o limiar de decisao volta a ser parametro explicito (``decision_threshold``),
    em vez do 0,5 fixo no corpo da funcao.

    B10: quando ``oof_per_repeat`` e fornecido, o IC da AUC passa a reamostrar
    conjuntamente PARTICIPANTES e REPETICOES da validacao cruzada. O IC da v23
    tratava as probabilidades medias como fixas e, por isso, subestimava a
    incerteza real, que inclui a variabilidade das particoes.
    """
    thr = float(acfg.decision_threshold)
    yhat = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yhat, labels=[0, 1]).ravel()
    auc = float(roc_auc_score(y, p))
    lo, hi = bootstrap_ci(lambda idx: float(roc_auc_score(y[idx], p[idx]))
                          if len(np.unique(y[idx])) == 2 else np.nan,
                          np.arange(len(y)), n_boot=acfg.n_boot_metrics,
                          seed=acfg.random_seed)
    out = {"n": int(len(y)), "n_mdd": int(y.sum()), "n_hc": int((1 - y).sum()),
           "decision_threshold": thr,
           "roc_auc": auc, "roc_auc_ci95_subjects_only": [lo, hi],
           "pr_auc": float(average_precision_score(y, p)),
           "balanced_accuracy": float(balanced_accuracy_score(y, yhat)),
           "sensitivity": float(tp / (tp + fn)) if (tp + fn) else float("nan"),
           "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
           "f1": float(f1_score(y, yhat)), "mcc": float(matthews_corrcoef(y, yhat)),
           "brier": float(brier_score_loss(y, p)),
           "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}}
    if oof_per_repeat is not None and np.ndim(oof_per_repeat) == 2:
        M = np.asarray(oof_per_repeat, float)
        n_rep, n_sub = M.shape
        rng = np.random.default_rng(acfg.random_seed)
        vals = []
        for _ in range(acfg.n_boot_metrics):
            ri = rng.integers(0, n_rep, n_rep)
            si = rng.integers(0, n_sub, n_sub)
            pb = np.nanmean(M[ri][:, si], axis=0)
            yb = y[si]
            if len(np.unique(yb)) == 2 and np.isfinite(pb).all():
                vals.append(roc_auc_score(yb, pb))
        if vals:
            out["roc_auc_ci95_subjects_and_cv"] = [float(np.percentile(vals, 2.5)),
                                                   float(np.percentile(vals, 97.5))]
            out["roc_auc_ci_note"] = ("IC que incorpora a variabilidade das particoes da CV; "
                                      "e o intervalo a reportar.")
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


def matched_observed_auc(X: np.ndarray, y: np.ndarray, acfg: AnalysisConfig) -> float:
    """B1: AUC observada sob o MESMO estimador usado na distribuicao nula."""
    r = nested_cv_scores(X, y, acfg, repeats=int(acfg.permutation_matched_repeats))
    return float(roc_auc_score(y, r["prob_mean"]))

def permutation_auc_test(X: np.ndarray, y: np.ndarray, observed_matched: float,
                         acfg: AnalysisConfig) -> Dict[str, Any]:
    """B1: teste de permutacao com ESTIMADOR IDENTICO entre observado e nulo.

    ERRO CORRIGIDO. Na v23 a AUC observada vinha de ``outer_repeats`` = 20
    repeticoes de CV (media de probabilidades OOF, portanto menos ruidosa) e cada
    AUC nula de UMA repeticao. Como a media sobre repeticoes altera a distribuicao
    amostral da estatistica, observado e nulo deixavam de ser permutaveis e o
    p-valor perdia interpretacao, com vies anticonservador.

    Agora ambos usam ``permutation_matched_repeats``. ``observed_matched`` DEVE ser
    calculado com o mesmo numero de repeticoes (ver ``matched_observed_auc``).
    """
    reps = int(acfg.permutation_matched_repeats)
    rng = np.random.default_rng(acfg.random_seed)
    null = []
    for _ in range(acfg.n_permutations):
        yp = rng.permutation(y)
        r = nested_cv_scores(X, yp, acfg, repeats=reps)
        null.append(float(roc_auc_score(yp, r["prob_mean"])))
    null = np.asarray(null, float)
    p = float((np.sum(null >= observed_matched) + 1) / (len(null) + 1))
    return {"n_permutations": int(len(null)),
            "estimator_repeats_observed_and_null": reps,
            "observed_auc_matched_estimator": float(observed_matched),
            "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
            "null_p95": float(np.percentile(null, 95)), "p_value": p,
            "note": ("Observado e nulo usam o mesmo numero de repeticoes de CV; "
                     "o p-valor NAO deve ser comparado com a AUC de %d repeticoes."
                     % acfg.outer_repeats)}



def environment_manifest() -> Dict[str, Any]:
    """Versoes de software que afetam a reproducibilidade numerica."""
    import scipy
    env = {"python": sys.version.split()[0], "platform": platform.platform(),
           "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
           "utc_created": datetime.now(timezone.utc).isoformat()}
    env["scikit_learn"] = __import__("sklearn").__version__ if _SKLEARN_OK else "not_installed"
    return env


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)



def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()



""" ==========================================
==============OBRIGATÓRIO MANTER==============
==============================================
"""
def fix_integer_wraparound(arr: np.ndarray, container_bits: int) -> Tuple[np.ndarray, int]:
    """Reinterpreta inteiros sem sinal como complemento de dois (v > 2**31 -> v - 2**32)."""
    half, full = 2 ** (container_bits - 1), 2 ** container_bits
    mask = arr > half
    n = int(mask.sum())
    if n == 0:
        return arr, 0
    out = arr.copy()
    out[mask] = out[mask] - full
    return out, n

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


def normalize_subject_id(values) -> pd.Series:
    """Normaliza identificadores para 8 digitos com zeros a esquerda."""
    s = pd.Series(values).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(8)



