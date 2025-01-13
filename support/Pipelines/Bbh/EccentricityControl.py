# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
from pathlib import Path
from typing import Optional, Sequence, Union

import click
import pandas as pd
import yaml

from spectre.Pipelines.Bbh.InitialData import generate_id
from spectre.Pipelines.EccentricityControl.EccentricityControlParams import (
    eccentricity_control_params,
    eccentricity_control_params_options,
)
from spectre.support.Schedule import scheduler_options

logger = logging.getLogger(__name__)


def eccentricity_control(
    h5_files: Union[Union[str, Path], Sequence[Union[str, Path]]],
    id_input_file_path: Union[str, Path],
    pipeline_dir: Union[str, Path],
    # Eccentricity control parameters
    tmin: Optional[float] = 500,
    tmax: Optional[float] = None,
    plot_output_dir: Optional[Union[str, Path]] = None,
    # Scheduler options
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

    - Get the new orbital parameters by calling the function
      'eccentricity_control_params' in
      'spectre.Pipelines.EccentricityControl.EccentricityControl'.

    - Displays the fit results in a tabular format using a pandas DataFrame.

    - If the eccentricity is below a threshold, it prints "Success" and
      indicates that the simulation can continue.

    - Generates new initial data based on updated orbital parameters using the
      'generate_id' function.

    Arguments:
      h5_files: files that contain the trajectory data
      id_input_file_path: path to the input file of the initial data run
      pipeline_dir : directory where the pipeline outputs are stored.

    See the 'eccentricity_control_params' function for details on the other
    arguments, as well as the 'schedule' function for the scheduling options.
    """
    # Read and process the initial data input file
    with open(id_input_file_path, "r") as open_input_file:
        id_metadata, id_input_file = yaml.safe_load_all(open_input_file)
    binary_data = id_input_file["Background"]["Binary"]
    id_params = id_metadata["Next"]["With"]
    control_params = id_params["control_params"]
    mass_A = control_params["mass_A"]
    mass_B = control_params["mass_B"]
    spin_A = control_params["spin_A"]
    spin_B = control_params["spin_B"]
    x_B, x_A = binary_data["XCoords"]
    separation = x_A - x_B

    # Find the current eccentricity and determine new parameters to put into
    # generate-id
    (
        eccentricity,
        ecc_std_dev,
        new_orbital_params,
    ) = eccentricity_control_params(
        h5_files,
        id_input_file_path,
        tmin=tmin,
        tmax=tmax,
        plot_output_dir=plot_output_dir,
    )

    # Create DataFrame to display data in tabular format
    data = {
        "Attribute": [
            "Eccentricity",
            "Eccentricity error",
            "Updated Omega0",
            "Updated adot0",
        ],
        "Value": [
            eccentricity,
            ecc_std_dev,
            new_orbital_params["Omega0"],
            new_orbital_params["adot0"],
        ],
    }
    df = pd.DataFrame(data)
    # Print header line
    print("=" * 40)
    # Display table
    print(df.to_string(index=False))
    print("=" * 40)

    # Stop eccentricity control if eccentricity is below threshold
    if eccentricity < 0.001:
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
        orbital_angular_velocity=new_orbital_params["Omega0"],
        radial_expansion_velocity=new_orbital_params["adot0"],
        # Scheduling options
        control=id_params["control"],
        refinement_level=id_params["control_refinement_level"],
        polynomial_order=id_params["control_polynomial_order"],
        evolve=True,
        eccentricity_control=True,
        pipeline_dir=pipeline_dir,
        **scheduler_kwargs,
    )


@click.command(name="eccentricity-control", help=eccentricity_control.__doc__)
@eccentricity_control_params_options
@click.option(
    "--pipeline-dir",
    "-d",
    type=click.Path(
        writable=True,
        path_type=Path,
    ),
    help="Directory where steps in the pipeline are created.",
)
@scheduler_options
def eccentricity_control_command(**kwargs):
    _rich_traceback_guard = True  # Hide traceback until here
    eccentricity_control(**kwargs)
