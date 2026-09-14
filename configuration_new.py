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

    electrode_labels = ['Fp1', 'Fpz', 'Fp2']




    # BLOCK BRAKE PARAMETERS
    block_time = 20.0 # em segundos


    # #FILTERS PARAMETERS    
    low_filter_frequency = 40.

    high_filter_frequency=1.

    band_limits = [high_filter_frequency,low_filter_frequency]

    filters_order = 4

    notch_quality_factor = 30 #Fator de Qualidade do Filtro Notch em dB

    notch_remove_frequency: Optional[float] = 50.0

    notch_mode: str = "auto"

    notch_line_ratio_threshold: float = 0.02

    # FEATURES PARAMETERS
    ocular_band: Tuple[float, float] = (0.5, 3.0)

    muscle_band: Tuple[float, float] = (20.0, 40.0)

    total_band: Tuple[float, float] = (1.0, 40.0)

    #Temporal PROTOCOL
    settling_ref_mode: str = "cohort"     # B3: cohort | self

    n_terminal_blocks: int = 5

    settle_tol_log2: float = 0.35

    settle_tol_js: float = 0.15

    tail_trend_warn: float = 0.50