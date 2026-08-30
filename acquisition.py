from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass(frozen=True)
class AcquisitionConfig:
    """Parametros de aquisicao do experimento pervasivo de 3 eletrodos."""
    fs: float = 250.0 # SAMPLING_FREQUENCY
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

    def __post_init__(self) -> None:
        if self.fs <= 0:
            raise ValueError("fs deve ser positivo.")
        if len(self.channel_names) != 3:
            raise ValueError("Esperados exatamente 3 canais.")
        if self.adc_bits > self.container_bits:
            raise ValueError("adc_bits nao pode exceder container_bits.")
