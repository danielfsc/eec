from dataclasses import dataclass
from typing import  Tuple, Optional
import numpy as np
@dataclass(frozen=True)
class Configuration:

    ## Acquisition Parameters
    
    sampling_frequency = 250 #Frequência da coleta de dados em Hz

    exclusion_time_window = 2.0

    bits_resolution = 24

    max_data_duration = 7200 # Duração do experimento em segundos

    number_electrodes = 3

    # #FILTERS PARAMETERS

    low_filter_frequency = 40.
    high_filter_frequency=1.
    band_limits = [high_filter_frequency,low_filter_frequency]
    filters_order = 4
    notch_quality_factor = 30 #Fator de Qualidade do Filtro Notch em dB
    notch_remove_frequency: Optional[float] = 50.0





    # #Schema Configurations
    # expected_n_channels: int = 3
    # require_finite: bool = True
    # require_integer: bool = True
    # max_duration_s: float = 7200.0
    # allowed_tasks: Tuple[str, ...] = ("still", "unspecified")
    # enforce_allowed_tasks: bool = True



    # #Acquisition Configurations
    # """Parametros de aquisicao do experimento pervasivo de 3 eletrodos."""
    # fs: float = 250.0
    # fs_is_assumption: bool = True   # B2: revertido para suposicao declarada
    # fs_verification_line_hz: Tuple[float, ...] = (50.0, 60.0)
    # fs_verification_min_prominence_log10: float = 0.30
    # fs_verification_tol_hz: float = 0.50
    # fs_source: str = ("Descritor MODMA: experimento pervasivo de 3 eletrodos, A/D de "
    #                       "24 bits a 250 Hz; lote unico verificado por esquema.")
    # protocol_id: str = "MODMA_3ch_pervasive_resting"
    # channel_names: Tuple[str, ...] = ("Fp1", "Fpz", "Fp2")
    # channel_order_source: str = "Ordem posicio" \
    # "nal das colunas do TXT (documentacao do dispositivo)."
    # container_bits: int = 32
    # adc_bits: int = 24
    # lsb_to_uv: Optional[float] = None

    # """Pre-processamento no sinal CONTINUO, antes de epocar."""
    # notch_hz: Optional[float] = 50.0
    # notch_q: float = 30.0
    # notch_mode: str = "auto"
    # notch_line_ratio_threshold: float = 0.02
    # bandpass_hz: Tuple[float, float] = (1.0, 40.0)
    # filter_order: int = 4
    # edge_trim_seconds: float = 2.0
    # epoch_seconds: float = 4.0
    # epoch_overlap: float = 0.0