from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class QCFreezeConfig:
    """B5: contrato de congelamento dos limiares de controle de qualidade.

    Limiares de QC nunca confrontados com a distribuicao empirica podem excluir
    zero ou quase toda a coorte sem que ninguem perceba antes de rodar. O
    procedimento obrigatorio e: (i) executar em modo integration; (ii) ler
    qc_distribution.csv, que informa quantos sujeitos CADA limiar excluiria;
    (iii) fixar os limiares; (iv) registrar aqui a data e a justificativa.
    Enquanto ``thresholds_frozen`` for False, o modo research fica bloqueado.
    """
    thresholds_frozen: bool = False
    freeze_date: Optional[str] = None
    freeze_source: str = ""
    require_freeze_for_research: bool = True

    def __post_init__(self) -> None:
        if self.thresholds_frozen and not self.freeze_date:
            raise ValueError("Limiares declarados congelados exigem freeze_date.")