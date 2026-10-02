from scipy.signal import butter,filtfilt, iirnotch
import numpy as np

from typing import Tuple

from configuration_new import Configuration

import utils2 as ut
default = Configuration()

def apply_filter(
        data:np.ndarray,
        cfg: Configuration  = default
        # sampling_frequency: float|int,
        # exclusion_time_window:float|int
    )->np.ndarray:

    data = subtract_mean_value(data)

    pre_ratio = line_noise_ratio(data, cfg.sampling_frequency, cfg.notch_remove_frequency or 50.0)

    if cfg.notch_mode == "always":
        do_notch = cfg.notch_remove_frequency is not None
    elif cfg.notch_mode == "never":
        do_notch = False
    else:
        do_notch = (cfg.notch_remove_frequency is not None) and (pre_ratio > cfg.notch_line_ratio_threshold)
    if do_notch:
        data = notch_filter(data)
    data = band_pass_filter(data)
    cut_size = int( round(cfg.sampling_frequency *cfg.exclusion_time_window))
    data = cut_array_edges(data, cut_size)
    return data , pos_filter_infos(data, do_notch, pre_ratio, cfg=cfg)

def band_pass_filter(data:np.ndarray,
                    band=default.band_limits,
                    sampling_frequency:int=default.sampling_frequency,
                    order=default.filters_order):

    nyq = 0.5 * sampling_frequency
    # print(band);

    numerator, denominator = butter(order, [ band[0]/nyq, band[1]/nyq], btype="band")

    return filtfilt(numerator, denominator, data, axis=-1)

def cut_array_edges(data:np.ndarray, size:int = 1):
    if(data.shape[1] <= 2*size ):
        raise ValueError("Sinal curto demais para o corte de borda.")
    return data[:, size:-size]

def fix_integer_wraparound(data: np.ndarray, bits_resolution: int = default.bits_resolution) -> Tuple[np.ndarray, int]:
    """
    Reinterpreta inteiros sem sinal como complemento de dois (v > 2**31 -> v - 2**32).
    """
    zero = 2 ** (bits_resolution - 1) 
    full = zero * 2
    mask = data > zero
    n = int(mask.sum())
    if n == 0:
        return data, 0
    data[mask] = data[mask] - full
    return data, n

def is_data_shape_ok(data:np.ndarray, cfg:Configuration=default):
    dur = data.shape[1] / cfg.sampling_frequency
    if dur > cfg.max_data_duration:
        raise ValueError("\n\tDado com duracao %.2f s acima do maximo.\n\t Verifique a variável \"max_data_duration\"" % ( dur))
    if data.shape[0] != cfg.number_electrodes:
        raise ValueError("\n\t Dado com número de eletrodos %d diferente do esperado.\n\t Verifique no arquivo de configuração a variável \"number_electrodes\"" % (data.shape[0]))
    return True

def line_noise_ratio(data: np.ndarray, sampling_frequency: float, notch_frequency: float = 50.0, bw: float = 1.0,
                     ref: Tuple[float, float] = (1.0, 45.0)) -> float:
    """Evidencia MEDIDA de interferencia de rede eletrica; base da decisao de notch."""
    f, p = ut.welch_power_spectral_density(data, sampling_frequency)
    return float(np.median(np.atleast_1d(ut.band_share(f, p, notch_frequency - bw, notch_frequency + bw, ref))))

def pos_filter_infos(data:np.ndarray, do_notch:bool, pre_ratio:float, cfg:Configuration = default):
    return {
        "notch_mode": cfg.notch_mode, 
        "notch_applied": bool(do_notch),
        "line_ratio_pre": float(pre_ratio),
        "line_ratio_post": float(line_noise_ratio(data, cfg.sampling_frequency, cfg.notch_remove_frequency or 50.0)),
        "edge_trim_s": float(cfg.exclusion_time_window),
        "n_samples_after_trim": int(data.shape[1])
    }

def subtract_mean_value(data:np.ndarray, axis=-1):
    return data - data.mean(axis=axis, keepdims=True)

def notch_filter(
        data:np.ndarray,
        remove_frequency = default.notch_remove_frequency,
        sampling_frequency = default.sampling_frequency, 
        quality_factor=default.notch_quality_factor ):
    
    numerator, denominator = iirnotch(remove_frequency, Q=quality_factor, fs=sampling_frequency)
    return  filtfilt(numerator, denominator, data, axis=-1)











