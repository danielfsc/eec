from dataclasses import dataclass
from typing import Optional, Tuple
@dataclass(frozen=True)
class PreprocConfig:
    """Pre-processamento no sinal CONTINUO, antes de epocar."""
    notch_hz: Optional[float] = 50.0
    notch_q: float = 30.0
    notch_mode: str = "auto"
    notch_line_ratio_threshold: float = 0.02
    bandpass_hz: Tuple[float, float] = (1.0, 40.0)
    filter_order: int = 4
    edge_trim_seconds: float = 2.0
    epoch_seconds: float = 4.0
    epoch_overlap: float = 0.0

    def __post_init__(self) -> None:
        lo, hi = self.bandpass_hz
        if not (0 < lo < hi):
            raise ValueError("bandpass invalida.")
        if self.notch_mode not in {"auto", "always", "never"}:
            raise ValueError("notch_mode deve ser auto|always|never.")
        if not (0.0 <= self.epoch_overlap < 1.0):
            raise ValueError("epoch_overlap fora de [0,1).")
