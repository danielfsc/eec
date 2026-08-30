from scipy.signal import butter,filtfilt, iirnotch
import numpy as np

from typing import Tuple
from configuration_new import Configuration

default = Configuration()

# def pass_filter(data:np.ndarray, cutoff, sampling_frequency, filter_type, filter_order, analog=False):
#     nyquist_frequency = 0.5 * sampling_frequency
#     normal_cutoff = cutoff / nyquist_frequency
#     b, a = butter(filter_order, normal_cutoff, btype=filter_type, analog=analog)
#     print([b,a])
#     filtered = filtfilt(b, a, data)
#     # print(data)
#     # print(y)
#     return filtered

# order = 4

# def high_pass_filter(data:np.ndarray, cutoff:float,sampling_frequency:int):
#     return pass_filter(data, cutoff , sampling_frequency, 'high', order)
    

# def low_pass_filter(data:np.ndarray, cutoff:float,sampling_frequency:int):
#     return pass_filter(data, cutoff , sampling_frequency, 'low', order)


def band_pass_filter(data:np.ndarray,
                    band=default.band_limits,
                    sampling_frequency:int=default.sampling_frequency,
                    order=default.filters_order):

    nyq = 0.5 * sampling_frequency
    # print(band);

    numerator, denominator = butter(order, [ band[0]/nyq, band[1]/nyq], btype="band")

    return filtfilt(numerator, denominator, data, axis=-1)


def notch_filter(
        data:np.ndarray,
        remove_frequency = default.notch_remove_frequency,
        sampling_frequency = default.sampling_frequency, 
        quality_factor=default.notch_quality_factor ):
    
    numerator, denominator = iirnotch(remove_frequency, Q=quality_factor, fs=sampling_frequency)
    return  filtfilt(numerator, denominator, data, axis=-1)


def subtract_mean_value(data:np.ndarray, axis=-1):
    return data - data.mean(axis=axis, keepdims=True)

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










