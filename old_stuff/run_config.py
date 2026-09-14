from dataclasses import dataclass, asdict, field
import hashlib
import json
from old_stuff.schema import SchemaConfig
from old_stuff.acquisition import AcquisitionConfig
from old_stuff.temporal_protocol import TemporalProtocolConfig
from old_stuff.pre_processor import PreprocConfig
from old_stuff.qc import QCConfig
from old_stuff.qc_freeze import QCFreezeConfig
from old_stuff.scale_inference import ScaleInferenceConfig
from old_stuff.cohort import CohortConfig
from old_stuff.feature import FeatureConfig
from old_stuff.analysis import AnalysisConfig

@dataclass(frozen=True)
class RunConfig:
    """Configuracao raiz imutavel com impressao digital reproduzivel."""
    acquisition: AcquisitionConfig = field(default_factory=AcquisitionConfig)
    schema: SchemaConfig = field(default_factory=SchemaConfig)
    temporal: TemporalProtocolConfig = field(default_factory=TemporalProtocolConfig)
    preproc: PreprocConfig = field(default_factory=PreprocConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    qc_freeze: QCFreezeConfig = field(default_factory=QCFreezeConfig)
    scale: ScaleInferenceConfig = field(default_factory=ScaleInferenceConfig)
    cohort: CohortConfig = field(default_factory=CohortConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    random_seed: int = 20240517

    def __post_init__(self) -> None:
        nyq = self.acquisition.fs / 2.0
        lo, hi = self.preproc.bandpass_hz
        if hi >= nyq:
            raise ValueError("Banda alta >= Nyquist.")
        if self.preproc.notch_hz is not None and self.preproc.notch_hz >= nyq:
            raise ValueError("notch acima de Nyquist.")

    def fingerprint(self) -> str:
        return hashlib.sha256(
            json.dumps(asdict(self), sort_keys=True, default=str).encode()).hexdigest()[:16]

