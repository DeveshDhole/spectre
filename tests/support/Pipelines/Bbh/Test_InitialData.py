# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
import shutil
import unittest
from pathlib import Path

import numpy.testing as npt
import yaml
from click.testing import CliRunner

from spectre.Informer import unit_test_build_path
from spectre.Pipelines.Bbh.InitialData import generate_id_command, id_parameters
from spectre.support.Logging import configure_logging


class TestInitialData(unittest.TestCase):
    def setUp(self):
        self.test_dir = Path(
            unit_test_build_path(), "support/Pipelines/Bbh/InitialData"
        )
        shutil.rmtree(self.test_dir, ignore_errors=True)
        self.test_dir.mkdir(parents=True, exist_ok=True)
        self.bin_dir = Path(unit_test_build_path(), "../../bin").resolve()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_generate_id(self):
        params = id_parameters(
            mass_a=0.6,
            mass_b=0.4,
            dimensionless_spin_a=[0.1, 0.2, 0.3],
            dimensionless_spin_b=[0.4, 0.5, 0.6],
            center_of_mass_offset=[0.1, 0.2, 0.3],
            linear_velocity=[0.1, 0.2, 0.3],
            separation=20.0,
            orbital_angular_velocity=0.01,
            radial_expansion_velocity=-1.0e-5,
            refinement_level=1,
            polynomial_order=5,
        )
        self.assertEqual(params["MassRight"], 0.6)
        self.assertEqual(params["MassLeft"], 0.4)
        self.assertEqual(params["XRight"], 8.0 + 0.1)
        self.assertEqual(params["XLeft"], -12.0 + 0.1)
        self.assertEqual(
            [params[f"CenterOfMassOffset_{yz}"] for yz in "yz"],
            [0.2, 0.3],
        )
        self.assertEqual(
            [params[f"LinearVelocity_{xyz}"] for xyz in "xyz"],
            [0.1, 0.2, 0.3],
        )
        self.assertAlmostEqual(params["ExcisionRadiusRight"], 1.07546791205)
        self.assertAlmostEqual(params["ExcisionRadiusLeft"], 0.5504049327)
        self.assertEqual(params["OrbitalAngularVelocity"], 0.01)
        self.assertEqual(params["RadialExpansionVelocity"], -1.0e-5)
        self.assertEqual(
            [params[f"DimensionlessSpinRight_{xyz}"] for xyz in "xyz"],
            [0.1, 0.2, 0.3],
        )
        self.assertEqual(
            [params[f"DimensionlessSpinLeft_{xyz}"] for xyz in "xyz"],
            [0.4, 0.5, 0.6],
        )
        npt.assert_allclose(
            [params[f"HorizonRotationRight_{xyz}"] for xyz in "xyz"],
            [-0.043236994315732, -0.086473988631464, -0.119710982947196],
        )
        npt.assert_allclose(
            [params[f"HorizonRotationLeft_{xyz}"] for xyz in "xyz"],
            [-0.337933017966707, -0.422416272458383, -0.49689952695006],
        )
        self.assertAlmostEqual(params["FalloffWidthRight"], 6.479672589667676)
        self.assertAlmostEqual(params["FalloffWidthLeft"], 5.520327410332324)
        self.assertEqual(params["L"], 1)
        self.assertEqual(params["P"], 5)
        # Newtonian center of mass (without offset) is zero
        self.assertAlmostEqual(
            params["MassRight"] * (params["XRight"] - 0.1)
            + params["MassLeft"] * (params["XLeft"] - 0.1),
            0.0,
        )

    def test_cli(self):
        common_args = [
            "--mass-ratio",
            "1.5",
            "--chi-A",
            "0.1",
            "0.2",
            "0.3",
            "--chi-B",
            "0.4",
            "0.5",
            "0.6",
            "--separation",
            "20",
            "--orbital-angular-velocity",
            "0.01",
            "--radial-expansion-velocity",
            "-1.0e-5",
            "--refinement-level",
            "1",
            "--polynomial-order",
            "5",
            "-E",
            str(self.bin_dir / "SolveXcts"),
            "--no-schedule",
        ]
        # Not using `CliRunner.invoke()` because it runs in an isolated
        # environment and doesn't work with MPI in the container.
        try:
            generate_id_command(
                common_args
                + [
                    "-o",
                    str(self.test_dir),
                    "--no-submit",
                ]
            )
        except SystemExit as e:
            self.assertEqual(e.code, 0)
        self.assertTrue(
            (self.test_dir / "ControlParams_000/InitialData.yaml").exists()
        )
        # Test with pipeline directory
        try:
            generate_id_command(
                common_args
                + [
                    "-d",
                    str(self.test_dir / "Pipeline"),
                    "--evolve",
                    "--eccentricity-control",
                    "--no-submit",
                ]
            )
        except SystemExit as e:
            self.assertEqual(e.code, 0)
        with open(
            self.test_dir
            / "Pipeline/000_InitialData/ControlParams_000/InitialData.yaml",
            "r",
        ) as open_input_file:
            metadata = next(yaml.safe_load_all(open_input_file))
        self.assertEqual(
            metadata["Next"],
            {
                "Run": "spectre.Pipelines.Bbh.PostprocessId:postprocess_id",
                "With": {
                    "id_input_file_path": "__file__",
                    "id_run_dir": "./",
                    "pipeline_dir": str(self.test_dir.resolve() / "Pipeline"),
                    "control": True,
                    "control_refinement_level": 1,
                    "control_polynomial_order": 5,
                    "control_params": {
                        "mass_A": 0.6,
                        "mass_B": 0.4,
                        "spin_A": [0.1, 0.2, 0.3],
                        "spin_B": [0.4, 0.5, 0.6],
                        "center_of_mass": [0.0, 0.0, 0.0],
                        "linear_momentum": [0.0, 0.0, 0.0],
                    },
                    "evolve": True,
                    "eccentricity_control": True,
                    "scheduler": "None",
                    "copy_executable": "None",
                    "submit_script_template": "None",
                    "submit": True,
                },
            },
        )


if __name__ == "__main__":
    configure_logging(log_level=logging.DEBUG)
    unittest.main(verbosity=2)
