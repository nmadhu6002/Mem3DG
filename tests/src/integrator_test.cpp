// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include <iostream>

#include <gtest/gtest.h>

#include "mem3dg/constants.h"
#include "mem3dg/mem3dg"
#include <Eigen/Core>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class IntegratorTest : public ::testing::Test {
public:
  IntegratorTest() {
    std::tie(mesh, vpg) = mem3dg::getIcosphereMatrix(1, 3);

    /// physical parameters
    p.bending.Kbc = 8.22e-5;
    p.tension.Ksg = 0.1;
    p.tension.At = 4.0 * mem3dg::constants::PI;
    p.osmotic.isPreferredVolume = true;
    p.osmotic.Kv = 0.01;
    p.osmotic.Vt = 4.0 / 3.0 * mem3dg::constants::PI * 0.7;
  }
  void SetUp() override {}
  void TearDown() override {}

  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> mesh;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vpg;

  mem3dg::solver::Parameters p;

  const double dt = 0.5, T = 50, eps = 0, tSave = 10;
  const std::string outputDir = "/tmp";
};

TEST_F(IntegratorTest, EulerIntegratorTest) {
  mem3dg::solver::System f(mesh, vpg, p, 0);
  f.initialize(0, true);
  mem3dg::solver::integrator::Euler integrator{f,   dt,        T, tSave,
                                               eps, outputDir, 0};
  integrator.trajFileName = "traj.nc";
  integrator.integrate();
}

TEST_F(IntegratorTest, ConjugateGradientIntegratorTest) {
  mem3dg::solver::System f(mesh, vpg, p, 0);
  f.initialize(0, true);
  mem3dg::solver::integrator::ConjugateGradient integrator{
      f, dt, T, tSave, eps, outputDir, 0};
  integrator.trajFileName = "traj.nc";
  integrator.integrate();
}

TEST_F(IntegratorTest, VelocityVerletIntegratorTest) {
  p.damping = 0.1 / dt;
  mem3dg::solver::System f(mesh, vpg, p, 0);
  f.initialize(0, true);
  mem3dg::solver::integrator::VelocityVerlet integrator{
      f, dt, 1, tSave, eps, outputDir, 0};
  integrator.trajFileName = "traj.nc";
  integrator.integrate();
}
