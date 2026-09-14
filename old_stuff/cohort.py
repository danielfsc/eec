from dataclasses import dataclass
from typing import Optional, Tuple
from vars import EXPECTED_SUBJECT_IDS_55

@dataclass(frozen=True)
class CohortConfig:
    """Definicao operacional da coorte-alvo e contagens tipo CONSORT."""
    expected_n: Optional[int] = 55
    expected_subject_ids: Optional[Tuple[str, ...]] = EXPECTED_SUBJECT_IDS_55
    duplicate_policy: str = "prefer_task"
    preferred_task: str = "still"
    require_one_file_per_subject: bool = True

