# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
import warnings
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import pandas as pd
import yaml

from spectre.Pipelines.Bbh.InitialData import generate_id
from spectre.Pipelines.EccentricityControl.EccentricityControl import (
    coordinate_separation_eccentricity_control,
)

# Suppress specific RuntimeWarnings
warnings.filterwarnings(
    "ignore", message="Number of calls to function has reached maxfev"
)

logger = logging.getLogger(__name__)


def eccentricity_control(
    h5_file,
    id_input_file_path,
    pipeline_dir: Union[str, Path],
    subfile_name_aha="ApparentHorizons/ControlSystemAhA_Centers.dat",
    subfile_name_ahb="ApparentHorizons/ControlSystemAhB_Centers.dat",
    tmin=500,
    tmax=None,  # use all available data
    output=None,  # reads from Inspiral.yaml
    **scheduler_kwargs,
):
    """Eccentricity reduction post inspiral.

    This function can be called after the inspiral has run (see the 'Next'
    section of the Inspiral.yaml file).

    This function does the following:

    - Reads orbital parameters from the 'id_input_file_path'.

    - Sets the time boundaries for the eccentricity reduction process, starting
      at 500 and using all available data by default, with the option to adjust
      'tmin' and 'tmax' dynamically.

    - Calls the 'coordinate_separation_eccentricity_control' function to
      calculate the current eccentricity and extract more accurate orbital
      parameters.

    - Displays the fit results in a tabular format using a pandas DataFrame.

    - If the eccentricity is below a threshold, it prints "Success" and
      indicates that the simulation can continue.

    - Generates new initial data based on updated orbital parameters using the
      'generate_id' function.

    Arguments:
      h5_file: file that contains the trajectory data
      id_input_file_path: path to the input file of the initial data run
      pipeline_dir : directory where the pipeline outputs are stored.
      subfile_name_aha: subfile for black hole A; optional parameter
      subfile_name_ahb: subfile for black hole B; optional parameter
      tmin: starting point for eccentricity reduction script; defaults to
        500 if not specified
      tmax: stopping point for eccentricity reduction script; set to use all
        available data by default
      output: outputs to terminal plus makes pdf file, if specified
    """
    # Read and process the initial data input file
    with open(id_input_file_path, "r") as open_input_file:
        id_metadata, id_input_file = yaml.safe_load_all(open_input_file)
    binary_data = id_input_file["Background"]["Binary"]
    orbital_angular_velocity = binary_data["AngularVelocity"]
    radial_expansion_velocity = binary_data["Expansion"]
    id_params = id_metadata["Next"]["With"]
    control_params = id_params["control_params"]
    mass_A = control_params["mass_A"]
    mass_B = control_params["mass_B"]
    spin_A = control_params["spin_A"]
    spin_B = control_params["spin_B"]
    x_B, x_A = binary_data["XCoords"]
    separation = x_A - x_B

    if output:
        fig = plt.figure()
    else:
        fig = None

    # Find the current eccentricity and extract new parameters to put into
    # generate-id. Call the coordinate_separation_eccentricity_control
    # function; extract only the ["H4"] , i.e. nonlinear fit.
    ecout = coordinate_separation_eccentricity_control(
        h5_file=h5_file,
        subfile_name_aha=subfile_name_aha,
        subfile_name_ahb=subfile_name_ahb,
        tmin=tmin,
        tmax=tmax,
        angular_velocity_from_xcts=orbital_angular_velocity,
        expansion_from_xcts=radial_expansion_velocity,
        fig=fig,
    )["H4"]

    if output:
        fig.savefig(output)

    # Print the fit result
    # extract new parameters to put into generate_id

    fit_result = ecout["fit result"]

    # Prepare data from fit result
    data = {
        "Attribute": [
            "Eccentricity",
            "Omega Update",
            "Expansion Update",
            "Previous Omega",
            "Previous Expansion",
            "Updated Omega",
            "Updated Expansion",
        ],
        "Value": [
            fit_result["eccentricity"],
            fit_result["xcts updates"]["omega update"],
            fit_result["xcts updates"]["expansion update"],
            fit_result["previous xcts values"]["omega"],
            fit_result["previous xcts values"]["expansion"],
            fit_result["updated xcts values"]["omega"],
            fit_result["updated xcts values"]["expansion"],
        ],
    }

    # Create DataFrame to display data in tabular format
    df = pd.DataFrame(data)
    # Print header line
    print("=" * 40)
    # Display table
    print(df.to_string(index=False))
    print("=" * 40)

    if fit_result["eccentricity"] < 0.001:
        print("Success")
        # Should continue the simulation either by restarting from a
        # checkpoint, or from the volume data - will do later
        return

    # Generate new initial data based on updated orbital parameters
    generate_id(
        mass_a=mass_A,
        mass_b=mass_B,
        dimensionless_spin_a=spin_A,
        dimensionless_spin_b=spin_B,
        # Orbital parameters
        separation=separation,
        orbital_angular_velocity=fit_result["updated xcts values"]["omega"],
        radial_expansion_velocity=fit_result["updated xcts values"][
            "expansion"
        ],
        # Scheduling options
        control=id_params["control"],
        refinement_level=id_params["control_refinement_level"],
        polynomial_order=id_params["control_polynomial_order"],
        evolve=True,
        eccentricity_control=True,
        pipeline_dir=pipeline_dir,
        **scheduler_kwargs,
    )
