// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <memory>
#include <pup.h>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Tags.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"

namespace {

void test_spherical_logical_coords() {
  for (size_t l = 2; l < 5; ++l) {
    const size_t nth = l + 1;
    const size_t nph = 2 * l + 1;
    const Mesh<2> mesh_s2{
        {nth, nph},
        {Spectral::Basis::SphericalHarmonic,
         Spectral::Basis::SphericalHarmonic},
        {Spectral::Quadrature::Gauss, Spectral::Quadrature::Equiangular}};
    const auto xi = logical_coordinates(mesh_s2);
    const ylm::Spherepack ylm(l, l);
    const auto xi_expected = ylm.theta_phi_points();
    CHECK(get<0>(xi) == xi_expected[0]);
    CHECK(get<1>(xi) == xi_expected[1]);
  }
}

SPECTRE_TEST_CASE("Unit.NumericalAlgorithms.Spectral.LogicalCoordinates",
                  "[Domain][Unit]") {
  test_spherical_logical_coords();
  using Affine2d =
      domain::CoordinateMaps::ProductOf2Maps<domain::CoordinateMaps::Affine,
                                             domain::CoordinateMaps::Affine>;
  using Affine3d =
      domain::CoordinateMaps::ProductOf3Maps<domain::CoordinateMaps::Affine,
                                             domain::CoordinateMaps::Affine,
                                             domain::CoordinateMaps::Affine>;

  const Mesh<1> mesh_1d{3, Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  const Mesh<2> mesh_2d{
      {{2, 3}}, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto};

  // [logical_coordinates_example]
  const Mesh<3> mesh_3d{{{5, 3, 2}},
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};

  const domain::CoordinateMaps::Affine x_map{-1.0, 1.0, -3.0, 7.0};
  const domain::CoordinateMaps::Affine y_map{-1.0, 1.0, -13.0, 47.0};
  const domain::CoordinateMaps::Affine z_map{-1.0, 1.0, -32.0, 74.0};

  const auto map_3d =
      domain::make_coordinate_map<Frame::ElementLogical, Frame::Grid>(
          Affine3d{x_map, y_map, z_map});

  const auto x_3d = map_3d(logical_coordinates(mesh_3d));
  // [logical_coordinates_example]

  const auto map_1d =
      domain::make_coordinate_map<Frame::ElementLogical, Frame::Grid>(
          domain::CoordinateMaps::Affine{x_map});
  const auto map_2d =
      domain::make_coordinate_map<Frame::ElementLogical, Frame::Grid>(
          Affine2d{x_map, y_map});
  const auto x_1d = map_1d(logical_coordinates(mesh_1d));
  const auto x_2d = map_2d(logical_coordinates(mesh_2d));

  CHECK(x_1d[0][0] == -3.0);
  CHECK(x_1d[0][1] == 2.0);
  CHECK(x_1d[0][2] == 7.0);

  CHECK(x_2d[0][0] == -3.0);
  CHECK(x_2d[0][1] == 7.0);

  CHECK(x_3d[0][0] == -3.0);
  CHECK(x_3d[0][2] == 2.0);
  CHECK(x_3d[0][4] == 7.0);

  CHECK(x_2d[1][0] == -13.0);
  CHECK(x_2d[1][2] == 17.0);
  CHECK(x_2d[1][4] == 47.0);

  CHECK(x_3d[1][0] == -13.0);
  CHECK(x_3d[1][5] == 17.0);
  CHECK(x_3d[1][10] == 47.0);

  CHECK(x_3d[2][0] == -32.0);
  CHECK(x_3d[2][15] == 74.0);

  TestHelpers::db::test_compute_tag<domain::Tags::LogicalCoordinates<1>>(
      "ElementLogicalCoordinates");
  TestHelpers::db::test_compute_tag<domain::Tags::LogicalCoordinates<2>>(
      "ElementLogicalCoordinates");
  TestHelpers::db::test_compute_tag<domain::Tags::LogicalCoordinates<3>>(
      "ElementLogicalCoordinates");
}
}  // namespace
