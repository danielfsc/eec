import pandas as pd
import numpy as np
import filter as ft
import utils2 as ut
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from configuration_new import Configuration

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

def data_load_filter_pre_analyze(
        data_path: str|Path,
        cfg: Configuration = Configuration()
    ) -> pd.DataFrame:
    """Carrega e filtra os dados de EEG de um arquivo CSV ou TXT.
    """
    # path_output =  f"output_{datetime.now():%Y_%m_%d_%H_%M_%S}/"
    # output = Path(path_output)
    # output_directories = {"integration":None, "quality":None, "analysis":None, "sensitivity":None}
    # for key in output_directories:
    #     output_directories[key]=output/key
    #     output_directories[key].mkdir(parents=True, exist_ok=True)
    metadata = ut.get_first_metadata(data_path)
    files, failures_files = ut.discover_eeg_txt_files(data_path)
    results = {}
    filtered_data = {}
    quality_control, profiles_rows, sampling_verification = [], [], []
    for file in files:
        try:
            subject_id, subject_state = file.name.split('.')[0].split('_')
            data_raw = pd.read_csv(file, sep=" ", names=cfg.electrode_labels)
            data_raw = np.transpose(data_raw.to_numpy().astype(float)) #Convertendo os dados para np.ndarray
            
            data_np, changed_values_on_wraparound = ft.fix_integer_wraparound(data_raw)
            data_np = np.ascontiguousarray(data_np)
            data_np, filter_info = ft.apply_filter(data_np, cfg) 

            filtered_data[subject_id] = data_np

            profile = ut.spectral_profile_blocked(data_np, cfg.sampling_frequency, cfg.block_time)
            profile[subject_id] = subject_id
            profile['task'] = subject_state
            profiles_rows.append(profile)

            metrics = ut.quality_control_metrics(data_raw, data_np, filtering_info=filter_info, subject_id=subject_id, file_name=file.name, eeg_task=subject_state, n_wraparound=changed_values_on_wraparound)

            quality_control.append(metrics)
            
        except Exception as error:
            print(f"Failed to process file {file.name}: {error}")
            failures_files.append({
                    "file": file.name, 
                    "subject_id": file.stem.split("_")[0],
                    "stage": "read_or_profile",
                    "error": "%s: %s" % (type(error).__name__, str(error)[:200])})

    return {
        "metadata": metadata,
        "profiles":profiles_rows,
        "quality_control": quality_control,
        "filtered_data": filtered_data,
        "failures_files": failures_files,
        "files": files
    }
