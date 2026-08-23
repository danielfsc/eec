
@dataclass(frozen=True)
class SubjectRecording:
    subject_id: str
    task: str
    data_counts: np.ndarray
    fs: float
    source_name: str
    parse_report: Dict[str, Any]
