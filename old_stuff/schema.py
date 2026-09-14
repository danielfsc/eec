from dataclasses import dataclass
from typing import  Tuple

@dataclass(frozen=True)
class SchemaConfig:
    """Contrato de esquema verificado em TODOS os arquivos antes de qualquer feature."""
    expected_n_channels: int = 3
    require_finite: bool = True
    require_integer: bool = True
    max_duration_s: float = 7200.0
    allowed_tasks: Tuple[str, ...] = ("still", "unspecified")
    enforce_allowed_tasks: bool = True

    