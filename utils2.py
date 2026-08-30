import numpy as np
import pandas as pd
from typing import Tuple

def convert_numpy_to_pandas(data:np.ndarray, columns:Tuple[str]|None):
    return pd.DataFrame(np.transpose(data), columns=columns)
