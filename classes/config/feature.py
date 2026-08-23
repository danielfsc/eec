
@dataclass(frozen=True)
class FeatureConfig:
    """Parametros PRE-ESPECIFICADOS de extracao (a janela vem de TemporalProtocolConfig)."""
    welch_nperseg_seconds: float = 2.0
    welch_overlap: float = 0.5
    bands: Tuple[Tuple[str, float, float], ...] = (
        ("delta", 1.0, 4.0), ("theta", 4.0, 8.0), ("alpha", 8.0, 13.0), ("beta", 13.0, 30.0))
    total_band: Tuple[float, float] = (1.0, 40.0)
    alpha_peak_band: Tuple[float, float] = (7.0, 13.0)
    alpha_min_peak_log10: float = 0.05
    aperiodic_fit_band: Tuple[float, float] = (3.0, 35.0)
    aperiodic_peak_sd: float = 1.5
    aperiodic_max_iter: int = 5
    ocular_band: Tuple[float, float] = (0.5, 3.0)
    muscle_band: Tuple[float, float] = (20.0, 40.0)
    perm_entropy_order: int = 3
    perm_entropy_delay: int = 1
    higuchi_kmax: int = 10
    aggregate: str = "median"
    min_good_epochs: int = 30

    def __post_init__(self) -> None:
        if self.aggregate not in ("median", "mean"):
            raise ValueError("aggregate deve ser median|mean.")
