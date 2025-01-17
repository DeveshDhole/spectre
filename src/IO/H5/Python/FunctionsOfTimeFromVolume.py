# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
from pathlib import Path
from typing import Optional, Union

import click
import numpy as np
import yaml
from rich.pretty import pretty_repr

import spectre.IO.H5 as spectre_h5
from spectre.Domain import deserialize_functions_of_time

logger = logging.getLogger(__name__)


def functions_of_time_from_volume(
    fot_vol_h5_path, fot_vol_subfile, match_time, fot_to_observe=None
):
    """This function returns a dictionary of the FunctionsOfTime from a volume
    h5 file evaluated at a desired time along with the time closest to the
    match_time passed in.

    Arguments:
    fot_vol_h5_path: The full path to volume data containing functions of time
    that have been evaluated at the match time.
    fot_vol_subfile: The subfile containing volume data with functions of time
    evaulated at match time.
    match_time: The desired time to retrieve the function of time values.
    """

    functions_of_time_at_match_time_dict = {}

    with spectre_h5.H5File(fot_vol_h5_path, "r") as h5file:
        if fot_vol_subfile.split(".")[-1] == "vol":
            fot_vol_subfile = fot_vol_subfile.split(".")[0]
        volfile = h5file.get_vol("/" + fot_vol_subfile)
        obs_ids = volfile.list_observation_ids()
        fot_times = np.array(list(map(volfile.get_observation_value, obs_ids)))
        which_obs_id = np.argmin(np.abs(fot_times - match_time))
        serialized_fots = volfile.get_functions_of_time(obs_ids[which_obs_id])
        functions_of_time = deserialize_functions_of_time(serialized_fots)
        logger.debug("Desired match time: " + str(match_time))
        logger.debug("Selected ObservationID: " + str(which_obs_id))
        logger.debug("Selected match time: " + str(fot_times[which_obs_id]))

        functions_of_time_at_match_time_dict["MatchTime"] = fot_times[
            which_obs_id
        ]

        for fot_name, fot in functions_of_time.items():
            fot_at_match_time = fot.func_and_2_derivs(fot_times[which_obs_id])
            if len(fot_at_match_time[0]) != 1:
                functions_of_time_at_match_time_dict[fot_name] = [
                    [coef for coef in x] for x in fot_at_match_time
                ]
            else:
                functions_of_time_at_match_time_dict[fot_name] = [
                    x[0] for x in fot_at_match_time
                ]

    return functions_of_time_at_match_time_dict
