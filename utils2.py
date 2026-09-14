import numpy as np
import pandas as pd
import scipy.signal as sg
from typing import Tuple, Dict, Any, List
from pathlib import Path

from vars import PROFILE_BANDS, NAME_REGEX
import filter as ft
from configuration_new import Configuration
default = Configuration()


def convert_numpy_to_pandas(data:np.ndarray, columns:Tuple[str]|None):
    return pd.DataFrame(np.transpose(data), columns=columns)



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


def rms_counts(data:np.ndarray):
    return float(np.median(np.sqrt((data ** 2).mean(axis=1))))


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



def welch_power_spectral_density(
        data: np.ndarray, 
        sampling_frequency: float|int, 
        time_window: float|int = 4.0):
    """PSD pelo método de Welch do sinal continuo (janelas de nperseg_s segundos)."""
    segment_length = min(int(round(time_window * sampling_frequency)), data.shape[-1])
    return sg.welch(data, fs=sampling_frequency, nperseg=segment_length, noverlap=segment_length // 2, detrend="constant", axis=-1)



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

def parse_modma_filename(name: str) -> Tuple[str, str]:
    """Extrai (subject_id, task) do nome do arquivo; task ausente -> 'unspecified'."""
    m = NAME_REGEX.match(name.strip())
    if not m:
        raise ValueError("Nome fora do padrão MODMA: %r" % name)
    return m.group("sid"), (m.group("task") or 'undefined')


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

