from dataclasses import dataclass

@dataclass(frozen=True)
class ScaleInferenceConfig:
    """Inferencia empirica de referencia e de escala (descritiva)."""
    enabled: bool = True
    anchor_rms_uv: float = 15.0
    anchor_source: str = "Valor central declarado a priori para RMS 1-40 Hz frontal de repouso."
    average_reference_tol: float = 1e-6
    common_mode_high: float = 0.50
    gcd_max_samples: int = 200000