
@dataclass(frozen=True)
class EvidencePolicy:
    """Contrato de evidencia. 'integration' nunca executa inferencia/ML."""
    mode: str = "integration"
    min_subjects_for_ml: int = 20
    min_per_class_for_ml: int = 20
    require_both_classes: bool = True
    forbid_clinical_scales: bool = True

    def __post_init__(self) -> None:
        if self.mode not in {"integration", "research"}:
            raise ValueError("mode deve ser 'integration' ou 'research'.")