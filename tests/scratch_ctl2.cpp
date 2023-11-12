#include <iostream>

#include <netcdf>

#include "mem3dg/mem3dg"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;

using EigenVectorX1d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
template <typename T>
using EigenVectorX1_T = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using EigenVectorX1i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3dr =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3ur =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3u = Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3>;
using EigenVectorX3sr =
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3s = Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>;

int main() {
    EigenVectorX3sr mesh;
    EigenVectorX3dr vpg;
    std::tie(mesh, vpg) = mem3dg::getIcosphereMatrix(1, 0);

    mem3dg::solver::Geometry geometry(mesh, vpg);

    mem3dg::solver::Parameters p;

    p.variation.isShapeVariation = true;

    p.temperature = 0;

    p.bending.Kb = 8.22e-5;
    p.bending.Kbc = 8.22e-5;
    p.bending.H0c = 10;



    mem3dg::solver::nProteinParameters p1;

    p1.isProteinVariation = true;
    p1.proteinMobility = 3;
    p1.epsilon = -1e-10;
    p1.eta = p.bending.Kb;

    std::vector<mem3dg::solver::nProteinParameters> pParameters = {p1};

    std::vector<EigenVectorX1d> pDensities = {Eigen::MatrixXd::Constant(vpg.rows(), 1, 0)};
    EigenVectorX3dr vel = Eigen::MatrixXd::Constant(vpg.rows(), 3, 0);

    mem3dg::solver::System sys(geometry, pDensities, vel, p, pParameters);

    const double dt = 0.25, T = 10, eps = 0, tSave = 10;
    std::size_t frame = 0;
    const std::string outputDir = "/tmp";

    mem3dg::solver::integrator::Euler integrator{sys, dt, T, tSave, eps, outputDir, frame};
    integrator.ifPrintToConsole = true;
    integrator.ifOutputTrajFile = true;
    integrator.integrate();

    return 0;
}
