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

#pragma once

#include <cassert>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <math.h>

#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

#include <vector>

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

struct Parameters {
  /// Bending modulus
  double Kb;
  /// Spontaneous curvature
  double H0;
  /// Sharpness of the spontaneous curvature hetergeneity
  double sharpness;
  /// radius of non-zero spontaneous curvature
  std::vector<double> r_H0;
  /// Global stretching modulus
  double Ksg;
  /// Vertex shifting constant
  double Kst;
  /// Local stretching modulus
  double Ksl;
  /// Edge spring constant
  double Kse;
  /// Volume regularization
  double Kv;
  /// Line tension
  double eta;
  /// binding energy per protein
  double epsilon;
  /// binding constant
  double Bc;
  /// Dissipation coefficient
  double gamma;
  /// Reduced volume
  double Vt;
  /// Ambient Pressure
  double cam;
  /// Temperature
  double temp;
  /// Noise
  double sigma;
  /// index of node with applied external force
  std::vector<double> pt;
  /// Magnitude of external force
  double Kf;
  /// level of concentration of the external force
  double conc;
  /// target height
  double height;
  /// domain of integration
  double radius;
  /// augmented Lagrangian parameter for area
  double lambdaSG = 0;
  /// augmented Lagrangian parameter for volume
  double lambdaV = 0;
};

struct Energy {
  /// total Energy of the system
  double totalE;
  /// kinetic energy of the membrane
  double kE;
  /// potential energy of the membrane
  double potE;
  /// bending energy of the membrane
  double BE;
  /// stretching energy of the membrane
  double sE;
  /// work of pressure within membrane
  double pE;
  /// chemical energy of the membrane protein
  double cE;
  /// line tension energy of interface
  double lE;
  /// work of external force
  double exE;
};

class DLL_PUBLIC System {
public:
  /// Parameters
  Parameters P;

  /// Cached mesh of interest
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh;
  gcs::ManifoldSurfaceMesh &mesh;
  /// Embedding and other geometric details
  std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
  gcs::VertexPositionGeometry &vpg;
  /// reference embedding geometry
  std::unique_ptr<gcs::VertexPositionGeometry> ptrrefVpg;
  gcs::VertexPositionGeometry &refVpg;
  /// Cached mesh data
  gcs::RichSurfaceMeshData richData;
  /// Energy
  Energy E;

  /// Cached bending stress
  gcs::VertexData<gc::Vector3> bendingPressure;
  /// Cached tension-induced capillary pressure
  gcs::VertexData<gc::Vector3> capillaryPressure;
  /// Cached interfacial line tension
  gcs::VertexData<gc::Vector3> lineTensionPressure;
  /// Cached externally-applied pressure
  gcs::VertexData<gc::Vector3> externalPressure;
  /// Cached relative inside pressure
  double insidePressure;
  /// Cached Surface tension
  double surfaceTension;

  /// Cached local stretching forces (in-plane regularization)
  gcs::VertexData<gc::Vector3> regularizationForce;
  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForce;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForce;

  /// Cached protein surface density
  gcs::VertexData<double> proteinDensity;
  /// Cached chemical potential
  gcs::VertexData<double> chemicalPotential;

  /// Whether or not do vertex shift
  const bool isVertexShift;
  /// Whether or not consider protein binding
  const bool isProtein;
  /// Whether adopt reduced volume parametrization
  const bool isReducedVolume;
  /// Whether calculate geodesic distance
  const bool isLocalCurvature;

  /// Target area per face
  gcs::FaceData<double> targetFaceAreas;
  /// Target total face area
  double targetSurfaceArea;
  /// Maximal volume
  double refVolume;
  /// Target length per edge
  gcs::EdgeData<double> targetEdgeLengths;
  /// Target edge cross length ratio
  gcs::EdgeData<double> targetLcr;
  /// Distance solver
  gcs::HeatMethodDistanceSolver heatSolver;

  /// Cached galerkin mass matrix
  Eigen::SparseMatrix<double> &M;
  /// Inverted galerkin mass matrix
  Eigen::SparseMatrix<double> M_inv;
  /// Cotangent Laplacian
  Eigen::SparseMatrix<double> &L;
  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromPtInd;

  /// L2 error norm
  double L2ErrorNorm;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// Interface Area;
  double interArea;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differencing past and current position
  gcs::VertexData<gc::Vector3> vel;
  // Mean curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H;
  // Gaussian curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> K;
  // Spontaneous curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H0;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// indices for vertices chosen for integration
  Eigen::Matrix<bool, Eigen::Dynamic, 1> mask;
  /// "the point" index
  size_t ptInd;

  /**
   * @brief Construct a new System object by reading unique_ptrs to mesh and
   * geometry object
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param richData_     Mesh rich data
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrrefVpg_, Parameters &p,
         bool isReducedVolume_, bool isProtein_, bool isLocalCurvature_,
         bool isVertexShift_) 

        :
        ptrmesh{ptrmesh_}, ptrvpg{ptrvpg_}, ptrrefVpg{ptrrefVpg_},
        mesh{}, vpg(), richData(), refVpg(), P(p),
        isReducedVolume(isReducedVolume_), isProtein(isProtein_),
        isLocalCurvature(isLocalCurvature_), isVertexShift(isVertexShift_),
        M(vpg.vertexLumpedMassMatrix), L(vpg.cotanLaplacian),
        bendingPressure(), insidePressure(0),
        capillaryPressure(),
        lineTensionPressure(), 
        chemicalPotential(),
        externalPressure(),
        regularizationForce(), targetLcr(),
        stochasticForce(), dampingForce(),
        proteinDensity(), vel(),
        E({0, 0, 0, 0, 0, 0, 0, 0, 0}), heatSolver()
         {
    
    mesh = *ptrmesh;
    vpg = *ptrvpg;
    refVpg = *ptrrefVpg;

    /// Initialize richData for ply file
    richData = gcs::RichSurfaceMeshData(*ptrMesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*ptrVpg);    
    richData = &richData;
    
    bendingPressure = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    capillaryPressure = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    lineTensionPressure = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    chemicalPotential = gcs::VertexData<double>(mesh, 0.0);
    externalPressure = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    regularizationForce = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    targetLcr = gcs::EdgeData<double>(mesh);
    stochasticForce = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    dampingForce = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    proteinDensity = gcs::VertexData<double>(mesh, 0);
    vel = gcs::VertexData<gc::Vector3>(mesh, {0, 0, 0});
    heatSolver(vpg);


    // GC computed properties
    vpg.requireFaceNormals();
    vpg.requireVertexLumpedMassMatrix();
    vpg.requireCotanLaplacian();
    vpg.requireFaceAreas();
    vpg.requireVertexIndices();
    vpg.requireVertexGaussianCurvatures();
    vpg.requireFaceIndices();
    vpg.requireEdgeLengths();
    vpg.requireVertexNormals();
    vpg.requireVertexDualAreas();
    vpg.requireCornerAngles();
    vpg.requireCornerScaledAngles();
    // vpg.requireVertexTangentBasis();

    // Check confliciting parameters and options
    checkParameters();

    // Initialize reference values
    initConstants();

    // Regularize the vetex position geometry if needed
    if (isVertexShift) {
      vertexShift(mesh, vpg, mask);
    }

    /// compute nonconstant values during simulation
    updateVertexPositions();
  }

  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   * Use case:
   * input: load two strings of mesh, one mesh is the initial configuration,
   * and the other the reference configuration, they the the same topology but
   * potentially different vertex position
   * output: tie of unique ptrs, both vpg associated with the same mesh
   *
   * @param inputMesh         Mesh connectivity
   * @param refMesh          Embedding and geometry information
   * @return tuple of mesh, vpg and refVpg
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshes(std::string inputMesh, std::string refMesh) {
    // assumes that the input and reference coordinates are using the same
    // mesh

    // Declare pointers to mesh / geometry objects
    std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
    std::unique_ptr<gcs::VertexPositionGeometry> inputCoordinates;
    std::unique_ptr<gcs::ManifoldSurfaceMesh> referenceMesh;
    std::unique_ptr<gcs::VertexPositionGeometry> referenceCoordinates;
    // gcs::VertexPositionGeometry *ptrRefVpg;
    size_t nSub = 0;

    // Load input mesh and geometry
    std::cout << "Loading input mesh " << inputMesh << " ...";
    std::tie(mesh, inputCoordinates) = gcs::readManifoldSurfaceMesh(inputMesh);
    std::cout << "Finished!" << std::endl;

    // Load input reference mesh and geometry
    std::cout << "Loading reference mesh " << refMesh << " ...";
    std::tie(referenceMesh, referenceCoordinates) =
        gcs::readManifoldSurfaceMesh(refMesh);
    std::cout << "Finished!" << std::endl;

    // Subdivide the mesh and geometry objects
    if (nSub > 0) {
      std::cout << "Subdivide input and reference mesh " << nSub
                << " time(s) ...";
      // mem3dg::subdivide(mesh, inputCoordinates, nSub);
      // mem3dg::subdivide(ptrRefMesh, ptrRefVpg, nSub);
      mem3dg::loopSubdivide(mesh, inputCoordinates, nSub);
      mem3dg::loopSubdivide(referenceMesh, referenceCoordinates, nSub);
      std::cout << "Finished!" << std::endl;
    }

    // Link referenceCoordinates to mesh instead of referenceMesh
    referenceCoordinates->mesh = *mesh;

    return std::make_tuple(std::move(mesh), std::move(inputCoordinates),
                           std::move(referenceCoordinates));
  }

  /**
   * @brief Construct a new System object by reading mesh path
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param richData_     Mesh rich data
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::string inputMesh, std::string refMesh, Parameters &p,
         bool isReducedVolume_, bool isProtein_, bool isLocalCurvature_,
         bool isVertexShift_)
      : System(readMeshes(inputMesh, refMesh), p, isReducedVolume_, isProtein_,
               isLocalCurvature_, isVertexShift_) {};

  /**
   * @brief Construct a new System object by reading tuple of unique_ptrs
   *
   * @param tuple        Mesh connectivity, Embedding and geometry
   * information, Mesh rich data
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local
   * curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             tuple,
         Parameters &p, bool isReducedVolume_, bool isProtein_,
         bool isLocalCurvature_, bool isVertexShift_)
      : System(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple), p,
               isReducedVolume_, isProtein_, isLocalCurvature_, isVertexShift_) {};

  // /**
  //  * @brief Construct a new System object by reference to mesh and geometry
  //  * object
  //  *
  //  * @param mesh_         Mesh connectivity
  //  * @param vpg_          Embedding and geometry information
  //  * @param richData_     Mesh rich data
  //  * @param p             Parameter of simulation
  //  * @param isReducedVolume Option of whether adopting reduced volume
  //  * parametrization
  //  * @param isProtein     Option of considering protein adsorption
  //  * @param isLocalCurvature Option of whether membrane has local curvature
  //  * @param isVertexShift Option of whether conducting vertex shift
  //  * regularization
  //  */
  // System(gcs::ManifoldSurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
  //        gcs::VertexPositionGeometry &refVpg_,
  //        gcs::RichSurfaceMeshData &richData_, Parameters &p,
  //        bool isReducedVolume_, bool isProtein_, bool isLocalCurvature_,
  //        bool isVertexShift_)
  //     : mesh(mesh_), vpg(vpg_), richData(richData_), refVpg(refVpg_), P(p),
  //       isReducedVolume(isReducedVolume_), isProtein(isProtein_),
  //       isLocalCurvature(isLocalCurvature_), isVertexShift(isVertexShift_),
  //       M(vpg.vertexLumpedMassMatrix), L(vpg.cotanLaplacian),
  //       bendingPressure(mesh_, {0, 0, 0}), insidePressure(0),
  //       capillaryPressure(mesh_, {0, 0, 0}),
  //       lineTensionPressure(mesh_, {0, 0, 0}), chemicalPotential(mesh_, 0.0),
  //       externalPressure(mesh_, {0, 0, 0}),
  //       regularizationForce(mesh_, {0, 0, 0}), targetLcr(mesh_),
  //       stochasticForce(mesh_, {0, 0, 0}), dampingForce(mesh_, {0, 0, 0}),
  //       proteinDensity(mesh_, 0), vel(mesh_, {0, 0, 0}),
  //       E({0, 0, 0, 0, 0, 0, 0, 0, 0}), heatSolver(vpg) {

  //   // GC computed properties
  //   vpg.requireFaceNormals();
  //   vpg.requireVertexLumpedMassMatrix();
  //   vpg.requireCotanLaplacian();
  //   vpg.requireFaceAreas();
  //   vpg.requireVertexIndices();
  //   vpg.requireVertexGaussianCurvatures();
  //   vpg.requireFaceIndices();
  //   vpg.requireEdgeLengths();
  //   vpg.requireVertexNormals();
  //   vpg.requireVertexDualAreas();
  //   vpg.requireCornerAngles();
  //   vpg.requireCornerScaledAngles();
  //   // vpg.requireVertexTangentBasis();

  //   // Check confliciting parameters and options
  //   checkParameters();

  //   // Initialize reference values
  //   initConstants();

  //   // Regularize the vetex position geometry if needed
  //   if (isVertexShift) {
  //     vertexShift(mesh, vpg, mask);
  //   }

  //   /// compute nonconstant values during simulation
  //   updateVertexPositions();
  // }

  /**
   * @brief Destroy the Force object
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~System() {
    vpg.unrequireFaceNormals();
    vpg.unrequireVertexLumpedMassMatrix();
    vpg.unrequireCotanLaplacian();
    vpg.unrequireFaceAreas();
    vpg.unrequireVertexIndices();
    vpg.unrequireVertexGaussianCurvatures();
    vpg.unrequireFaceIndices();
    vpg.unrequireEdgeLengths();
    vpg.unrequireVertexNormals();
    vpg.unrequireVertexDualAreas();
    vpg.unrequireCornerAngles();
    vpg.unrequireCornerScaledAngles();
  }

  // ==========================================================
  // ================     Initialization     ==================
  // ==========================================================
  /**
   * @brief Check all conflicting parameters and options
   *
   */
  void checkParameters();

  /**
   * @brief testing of random number generator pcg
   *
   */
  void pcg_test();

  /**
   * @brief Initialize all constant values (on refVpg) needed for computation
   *
   */
  void initConstants();

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   *
   */
  void updateVertexPositions();

  // ==========================================================
  // ================        Pressure        ==================
  // ==========================================================
  /**
   * @brief Get bending pressure component of the system
   */
  void getBendingPressure();

  /**
   * @brief Get chemical potential of the system
   */
  void getChemicalPotential();

  /**
   * @brief Get capillary pressure component of the system
   */
  void getCapillaryPressure();

  /**
   * @brief Get inside pressure component of the system
   */
  void getInsidePressure();

  /**
   * @brief Get regularization pressure component of the system
   */
  void getRegularizationForce();

  /**
   * @brief Get line tension pressure component of the system
   */
  void getLineTensionPressure();

  /**
   * @brief Get DPD forces of the system
   */
  void getDPDForces();

  /**
   * @brief Get external pressure component of the system
   */
  void getExternalPressure();

  /**
   * @brief Get all forces of the system
   */
  void getAllForces();

  // ==========================================================
  // ================        Energy          ==================
  // ==========================================================
  /**
   * @brief Get bending energy
   */
  void getBendingEnergy();

  /**
   * @brief Get surface energy
   */
  void getSurfaceEnergy();

  /**
   * @brief Get pressure work
   */
  void getPressureEnergy();

  /**
   * @brief Get chemical energy
   */
  void getChemicalEnergy();

  /**
   * @brief Get line tension energy
   */
  void getLineTensionEnergy();

  /**
   * @brief Get external force energy
   */
  void getExternalForceEnergy();

  /**
   * @brief Get kinetic energy
   */
  void getKineticEnergy();

  /**
   * @brief Get potential energy
   */
  void getPotentialEnergy();

  /**
   * @brief Get all components of energy (free energy)
   */
  void getFreeEnergy();

  /**
   * @brief Get the L2 Error norm of the PDE
   */
  void getL2ErrorNorm(Eigen::Matrix<double, Eigen::Dynamic, 3> pressure);

  /**
   * @brief Get the L2 norm of the pressure
   */
  double getL2Norm(Eigen::Matrix<double, Eigen::Dynamic, 3> pressure) const;
};
} // namespace mem3dg
