# Function to get events from snirf file
from collections import Counter
from pathlib import Path
from typing import Union, Optional

import h5py
import numpy as np
import pandas as pd

COL_TIMESTAMP = "Timestamp"
COL_DURATION = "Duration"
COL_VALUE = "Value"
COL_EVENT = "Event"

COLS_DEFAULT = [COL_TIMESTAMP, COL_DURATION, COL_VALUE]


def get_events_from_snirf(filename: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Read events from a SNIRF file as a dataframe

    Parameters
    ----------
    filename: Union[str, Path]
        the path to the SNIRF file

    Returns
    -------
    pd.DataFrame
        a pandas dataframe with events, sorted by timestamp.
    """
    df_events = pd.DataFrame()

    with h5py.File(Path(filename).resolve(), "r") as file:
        event_dfs = [
            pd.DataFrame(
                np.array(stim_info["data"]),
                columns=(np.array(stim_info["dataLabels"]).astype(str) if "dataLabels" in stim_info else COLS_DEFAULT),
            ).assign(Event=[np.array(stim_info["name"]).astype(str)] * len(np.array(stim_info["data"])))
            for stim_name, stim_info in sorted(file["nirs"].items())
            if stim_name.startswith("stim")
        ]

    if not event_dfs:
        return
    df_events = pd.concat(event_dfs, ignore_index=True)

    df_events.sort_values(by=[COL_TIMESTAMP], inplace=True, ignore_index=True)
    df_events.reset_index(inplace=True, drop=True)
    df_events = df_events[
        [COL_TIMESTAMP, COL_EVENT, COL_DURATION]
        + [col for col in df_events.columns if col not in [COL_TIMESTAMP, COL_EVENT, COL_DURATION]]
    ]

    # merge one-hot encoded columns
    col_prefixes = [col.split(".")[0] for col in df_events.columns]
    counter = Counter(col_prefixes)

    for col_prefix in set(col_prefixes):
        if counter[col_prefix] == 1:
            continue
        cols = [col for col in df_events.columns if col.startswith(col_prefix)]
        for col in cols:
            df_events.loc[df_events[col] == 1.0, col_prefix] = ".".join(col.split(".")[1:])
        df_events.drop(columns=cols, inplace=True)

    return df_events

