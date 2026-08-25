from dataclasses import dataclass
from typing import  Tuple, Optional

@dataclass(frozen=True)
class Configuration:
    #Schema Configurations
    expected_n_channels: int = 3
    require_finite: bool = True
    require_integer: bool = True
    max_duration_s: float = 7200.0
    allowed_tasks: Tuple[str, ...] = ("still", "unspecified")
    enforce_allowed_tasks: bool = True



    #Acquisition Configurations
    """Parametros de aquisicao do experimento pervasivo de 3 eletrodos."""
    fs: float = 250.0
    fs_is_assumption: bool = True   # B2: revertido para suposicao declarada
    fs_verification_line_hz: Tuple[float, ...] = (50.0, 60.0)
    fs_verification_min_prominence_log10: float = 0.30
    fs_verification_tol_hz: float = 0.50
    fs_source: str = ("Descritor MODMA: experimento pervasivo de 3 eletrodos, A/D de "
                          "24 bits a 250 Hz; lote unico verificado por esquema.")
    protocol_id: str = "MODMA_3ch_pervasive_resting"
    channel_names: Tuple[str, ...] = ("Fp1", "Fpz", "Fp2")
    channel_order_source: str = "Ordem posicio" \
    "nal das colunas do TXT (documentacao do dispositivo)."
    container_bits: int = 32
    adc_bits: int = 24
    lsb_to_uv: Optional[float] = None