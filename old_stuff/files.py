from pathlib import Path
from typing import Tuple, List, Dict, Any
import pandas as pd
import numpy as np
import warnings

from vars import NAME_REGEX, TASK_UNSPECIFIED
from old_stuff.configuration import Configuration
from old_stuff.subject_recording import SubjectRecording
import old_stuff.utils as ut

def list_data_eeg_files(data_dir: str | Path) -> Tuple[List[Path], List[Dict[str, str]]]:
    """Retorna os arquivos em um diretório:
    INPUT:
    data_dir: str|Path - caminho para o diretório com entrada de dados

    OUTPUT:
    valid_name: List[Path] - Arquivos com o nome no padrão modma %sid%_%task%.txt
    invalid_name:List[Path] => Arquivos que estão com o nome fora do padrão.
    """
    root = Path(data_dir)
    if not root.is_dir():
        raise FileNotFoundError("Diretório com dados EEG inexistente: %s" % root)
    valid_name, invalid_name = [], []
    for p in sorted(root.glob("*.txt")):
        try:
            get_modma_sid_task(p.name)
            valid_name.append(p)
        except ValueError as exc:
            invalid_name.append({"file": p.name, "stage": "filename", "error": str(exc)[:200]})
    if not valid_name:
        raise FileNotFoundError("Nenhum TXT MODMA valido encontrado.")
    return valid_name, invalid_name

def get_modma_sid_task(name: str) -> Tuple[str, str]:
    """Extrai (subject_id, task) do nome do arquivo; task ausente -> 'unspecified'.
    """
    m = NAME_REGEX.match(name.strip())
    if not m:
        raise ValueError("Nome fora do padrão MODMA: %r" % name)
    return m.group("sid"), (m.group("task") or TASK_UNSPECIFIED)


def load_modma_metadata(path: str | Path) -> pd.DataFrame:
    """Le a planilha, canoniza IDs, mapeia MDD/HC e extrai covariáveis demográficas e retorna um dataframe do pandas
    INPUT:
    path: str|Path - para o arquivo de metadados. O arquivo tem que estar em XLS, XLSX ou CSV
    OUTPUT: Pandas.Dataframe
    """
    p = Path(path)
    raw = pd.read_csv(p) if p.suffix.lower() == ".csv" else pd.read_excel(p)
    if raw.empty:
        raise ValueError("Planilha de metadados vazia.")
    raw = raw.loc[:, [c for c in raw.columns if not str(c).lower().startswith("unnamed")]]
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]
    sidcols = [c for c in raw.columns
               if c in {"subject_id", "subject", "id", "participant_id"} or c.startswith("subject")]
    groupcols = [c for c in raw.columns
                 if c in {"type", "group", "class", "diagnosis", "diagnostic_group"}]
    if len(sidcols) != 1 or len(groupcols) != 1:
        raise ValueError("Esperada exatamente uma coluna de ID e uma de grupo.")
    out = raw.rename(columns={sidcols[0]: "subject_id", groupcols[0]: "label_raw"}).copy()
    out["subject_id"] = ut.normalize_subject_id(out["subject_id"])
    if out["subject_id"].duplicated().any():
        raise ValueError("IDs duplicados nos metadados.")
    lab = out["label_raw"].astype(str).str.upper().str.strip()
    out["label"] = np.where(lab.str.startswith("MDD"), 1, np.where(lab.str.startswith("HC"), 0, -1))
    if (out["label"] < 0).any():
        raise ValueError("Rotulo nao mapeável para MDD/HC.")
    for c in list(out.columns):
        if c.startswith("education"):
            out = out.rename(columns={c: "education_years"})
        elif c.startswith("gender") or c.startswith("sex"):
            out["sex"], sex_report = ut.parse_sex_column(out[c])
    if "sex" in out.columns and not out["sex"].notna().any():
        warnings.warn("B9: coluna de sexo presente mas totalmente nao mapeada; os modelos "
                      "ajustados por confundidores ficariam silenciosamente vazios.",
                      RuntimeWarning)
    out["id_prefix"] = out["subject_id"].str[:4]
    keep = [c for c in ("subject_id", "label", "label_raw", "age", "sex",
                        "education_years", "id_prefix") if c in out.columns]
    return out[keep].drop_duplicates("subject_id").reset_index(drop=True)

# def load_modma_file(path: str | Path, acq: AcquisitionConfig, source_name: str, schema: SchemaConfig) -> SubjectRecording:

def load_modma_file(path: str | Path, cfg: Configuration, source_name: str) -> SubjectRecording:
    """Le um TXT MODMA de 3 canais e valida o esquema antes de devolver o registro."""
    sid, task = get_modma_sid_task(source_name)
    arr = np.loadtxt(str(path), dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != cfg.expected_n_channels:
        raise ValueError("%s: esperado (N,%d), obtido %s."
                         % (source_name, cfg.expected_n_channels, arr.shape))
    if cfg.require_finite and not np.isfinite(arr).all():
        raise ValueError("%s: valores nao finitos." % source_name)
    if cfg.require_integer and not np.allclose(arr, np.round(arr)):
        raise ValueError("%s: valores nao inteiros." % source_name)
    if cfg.enforce_allowed_tasks and task not in cfg.allowed_tasks:
        raise ValueError("%s: tarefa '%s' fora do protocolo." % (source_name, task))
    arr, n_wrap = ut.fix_integer_wraparound(arr, cfg.container_bits)
    x = np.ascontiguousarray(arr.T)
    dur = x.shape[1] / cfg.fs
    if dur > cfg.max_duration_s:
        raise ValueError("%s: duracao %.1f s acima do maximo." % (source_name, dur))
    rep = {"n_samples": int(x.shape[1]), "duration_s": float(dur),
           "n_wraparound_fixed": int(n_wrap),
           "wraparound_fraction": float(n_wrap / x.size),
           "dc_offset_counts": [float(v) for v in x.mean(axis=1)],
           "task_parsed": task}
    return SubjectRecording(sid, task, x, cfg.fs, source_name, rep)

