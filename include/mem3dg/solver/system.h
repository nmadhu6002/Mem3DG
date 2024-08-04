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

// #include <cassert>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <functional>
#include <iomanip>
#include <math.h>
#include <vector>

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/type_utilities.h"

#include "mem3dg/solver/forces.h"
#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/parameters.h"
#include "mem3dg/solver/nparameters.h"
#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/mutable_trajfile.h"
#include "mem3dg/solver/trajfile.h"
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/// @brief mem3dg namespace
namespace mem3dg {
/// @brief solver namespace
namespace solver {

struct Energy {
  /// time
  double time = 0;
  /// total Energy of the system
  double totalEnergy = 0;
  /// kinetic energy of the membrane
  double kineticEnergy = 0;
  /// potential energy of the membrane
  double potentialEnergy = 0;
  /// spontaneous curvature energy of the membrane
  double spontaneousCurvatureEnergy = 0;
  /// Gaussian curvature energy
  double gaussianCurvatureEnergy = 0;
  /// deviatoric curvature energy of the membrane
  double deviatoricCurvatureEnergy = 0;
  /// area difference energy of the membrane
  double areaDifferenceEnergy = 0;
  /// stretching energy of the membrane
  double surfaceEnergy = 0;
  /// work of pressure within membrane
  double pressureEnergy = 0;
  /// adsorption energy of the membrane protein
  double adsorptionEnergy = 0;
  /// aggregation energy of the membrane protein
  double aggregationEnergy = 0;
  /// entropy energy of the membrane protein
  double entropyEnergy = 0;
  /// line tension energy of interface
  double dirichletEnergy = 0;
  /// work of external force
  double externalWork = 0;
  /// protein interior penalty energy
  double proteinInteriorPenalty = 0;
  /// membrane self-avoidance penalty energy
  double selfAvoidancePenalty = 0;
  /// mesh edge spring energy
  double edgeSpringEnergy = 0;
  /// mesh face spring energy
  double faceSpringEnergy = 0;
  /// mesh LCR spring energy
  double lcrSpringEnergy = 0;
};

class DLL_PUBLIC System {
protected:
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;

public:
  /// Parameters
  Parameters parameters;
  /// Multiple Protein Parameters
  std::vector<nProteinParameters> pParameters;
  /// Mesh processor
  MeshProcessor meshProcessor;
  /// Geometry
  Geometry &geometry;
  /// Energy
  Energy energy;
  /// Time
  double time;
  /// Forces of the system
  Forces forces;

  /// mechanical error norm
  double mechErrorNorm;
  /// chemical error norm
  double chemErrorNorm;
  std::vector<double> chemErrorNorms;
  /// Cached protein surface density
  gcs::VertexData<double> proteinDensity;
  /// Spontaneous curvature gradient of the mesh
  gcs::FaceData<gc::Vector3> proteinDensityGradient;

  std::vector<gcs::VertexData<double>> pDensities;
  std::vector<gcs::FaceData<gc::Vector3>> pDensityGradients;
  std::vector<std::vector<double>> previousAreas;

  /// Cached vertex velocity
  gcs::VertexData<gc::Vector3> velocity;
  /// Cached vertex protein rate of change
  gcs::VertexData<double> proteinRateOfChange;

  std::vector<gcs::VertexData<double>> pRatesOfChange;

  /// Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Bending rigidity of the membrane
  gcs::VertexData<double> Kb;
  /// deviatoric rigidity of the membrane
  gcs::VertexData<double> Kd;
  /// gaussian modulus of the membrane
  gcs::VertexData<double> Kg;
  /// is Smooth
  bool isSmooth;
  /// if being mutated
  gcs::VertexData<bool> mutationMarker;
  /// projected time of collision
  double projectedCollideTime;

  // =======================================
  // =======       NetCDF Files     ========
  // =======================================
#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param trajFile NetCDF trajectory file
   * @param p Parameters of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_, std::string trajFile, int startingFrame,
         Parameters &p)
      : System(geometry_, readTrajFile(trajFile, startingFrame), p){};
  // System(std::string trajFile, int startingFrame, Parameters &p)
  //     : System(readTrajFile(trajFile, startingFrame), p){};
  System(Geometry &geometry_, std::string trajFile, int startingFrame,
         Parameters &p, std::vector<nProteinParameters> pParameters_)
      : System(geometry_, nProteinReadTrajFile(trajFile, startingFrame, pParameters_.size()), p, pParameters_){};
#endif

private:
  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param tuple <vertex protein density, velocity of the vertex, time of the
   * system>
   * @param p Parameters of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_,
         std::tuple<EigenVectorX1d, EigenVectorX3dr, double> tuple,
         Parameters &p)
      : System(geometry_, std::get<0>(tuple), std::get<1>(tuple), p,
               std::get<2>(tuple)) {}

  System(Geometry &geometry_,
         std::tuple<std::vector<EigenVectorX1d>, EigenVectorX3dr, double> tuple,
         Parameters &p, std::vector<nProteinParameters> pParameters_)
      : System(geometry_, std::get<0>(tuple), std::get<1>(tuple), p, pParameters_,
               std::get<2>(tuple)) {}

public:
  // =======================================
  // =======    Geometry Central    ========
  // =======================================

  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param proteinDensity_ vertex protein density
   * @param velocity_ velocity of the vertex
   * @param p parameters struct of the system
   * @param time_ time of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_, EigenVectorX1d &proteinDensity_,
         EigenVectorX3dr &velocity_, Parameters &p, double time_ = 0)
      : System(geometry_, proteinDensity_, velocity_, time_) {
    parameters = p;
  }

  System(Geometry &geometry_, std::vector<EigenVectorX1d> pDensities_,
         EigenVectorX3dr &velocity_, Parameters &p,
         std::vector<nProteinParameters> pParameters_, double time_ = 0)
      : geometry(geometry_), forces(geometry, pDensities_.size()), time(time_) {
    int n = pDensities_.size();
    energy = Energy({time, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    proteinDensity = gc::VertexData<double>(*geometry.mesh, 0);
    proteinDensityGradient =
        gcs::FaceData<gc::Vector3>(*geometry.mesh, {0, 0, 0});
    velocity = gcs::VertexData<gc::Vector3>(*geometry.mesh, {0, 0, 0});
    proteinRateOfChange = gcs::VertexData<double>(*geometry.mesh, 0);
    H0 = gcs::VertexData<double>(*geometry.mesh);
    Kb = gcs::VertexData<double>(*geometry.mesh);
    Kd = gcs::VertexData<double>(*geometry.mesh);
    Kg = gcs::VertexData<double>(*geometry.mesh);

    chemErrorNorm = 0;
    mechErrorNorm = 0;

    isSmooth = true;
    mutationMarker = gc::VertexData<bool>(*geometry.mesh, false);
    parameters = p;

    for (int i = 0; i < n; ++i) {
      chemErrorNorms.push_back(0);

      pDensities.push_back(gc::VertexData<double>(*geometry.mesh, 1));
      pDensityGradients.push_back(gcs::FaceData<gc::Vector3>(*geometry.mesh, {0, 0, 0}));
      pRatesOfChange.push_back(gcs::VertexData<double>(*geometry.mesh, 0));

      pDensities[i].raw() = pDensities_[i];

      nProteinParameters blank;
      pParameters.push_back(blank);

      pParameters[i] = pParameters_[i];
    }
    toMatrix(velocity) = velocity_;
  }

  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param proteinDensity_ vertex protein density
   * @param velocity_ velocity of the vertex
   * @param time_ time of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_, EigenVectorX1d &proteinDensity_,
         EigenVectorX3dr &velocity_, double time_ = 0)
      : System(geometry_, time_) {
    proteinDensity.raw() = proteinDensity_;
    toMatrix(velocity) = velocity_;
  }

  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param p parameters struct of the system
   * @param time_ time of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_, Parameters &p, double time_ = 0)
      : System(geometry_, time_) {
    parameters = p;
  }

  /**
   * @brief Construct System
   * @param geometry_ geometry of the system
   * @param time_ time of the system
   * @return system instance
   *
   */
  System(Geometry &geometry_, double time_ = 0)
      : geometry(geometry_), time(time_), forces(geometry) {
    energy = Energy({time, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    proteinDensity = gc::VertexData<double>(*geometry.mesh, 1);
    proteinDensityGradient =
        gcs::FaceData<gc::Vector3>(*geometry.mesh, {0, 0, 0});
    velocity = gcs::VertexData<gc::Vector3>(*geometry.mesh, {0, 0, 0});
    proteinRateOfChange = gcs::VertexData<double>(*geometry.mesh, 0);
    H0 = gcs::VertexData<double>(*geometry.mesh);
    Kb = gcs::VertexData<double>(*geometry.mesh);
    Kd = gcs::VertexData<double>(*geometry.mesh);
    Kg = gcs::VertexData<double>(*geometry.mesh);

    chemErrorNorm = 0;
    mechErrorNorm = 0;

    isSmooth = true;
    mutationMarker = gc::VertexData<bool>(*geometry.mesh, false);
  }

public:
  /**
   * @brief Destroy the System
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~System() {}

  // ==========================================================
  // ================          io.cpp        ==================
  // ==========================================================

  /**
   * @brief Map the continuation variables
   *
   */
  void mapContinuationVariables(std::string plyFile);

  /**
   * @brief Save RichData to .ply file
   *
   */
  void saveRichData(std::string PathToSave, bool isJustGeometry = false);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct a tuple of unique_ptrs from netcdf path
   *
   */
  // std::tuple<Geometry &&, EigenVectorX1d &, EigenVectorX3dr &, double>
  // readTrajFile(std::string trajFile, int startingFrame);
  std::tuple<EigenVectorX1d, EigenVectorX3dr, double>
  readTrajFile(std::string trajFile, int startingFrame);
  std::tuple<std::vector<EigenVectorX1d>, EigenVectorX3dr, double>
  nProteinReadTrajFile(std::string trajFile, int startingFrame, int n);
#endif

  // ==========================================================
  // ================     init.cpp           ==================
  // ==========================================================
  /**
   * @brief Check all conflicting parameters and options
   *
   */
  void checkConfiguration();

  /**
   * @brief Initialize system
   *
   */
  void initialize(bool ifMutateMesh = false);

  /**
   * @brief Update the vertex position and recompute cached values
   *
   * All vpg dependent quantities, volume, and area of mesh
   *
   * @warning Function may invalidate the memory address of existing Eigen
   * Mapped views!
   */
  void updateConfigurations();

  /**
   * @brief update various prescription of based on scalar and functional
   * parameters
   */
  bool updatePrescription(std::map<std::string, double> &lastUpdateTime,
                          double timeStep);
  bool
  updatePrescription(bool &ifMutateMesh, bool &ifUpdateNotableVertex,
                     std::vector<bool> &ifUpdateNotableVertices, bool &ifUpdateGeodesics,
                     bool &ifUpdateProteinDensityDistribution,
                     std::vector<bool> &ifUpdateProteinDensityDistributions,
                     bool &ifUpdateMask);

  // ==========================================================
  // ================        Force.cpp       ==================
  // ==========================================================
  /**
   * @brief Compute and update all conservative forces, update
   * mechanicalForce(Vec) with conservativeForce(Vec)
   */
  void computeConservativeForcing();

  /**
   * @brief Compute and append all non-conservative forces, update
   * mechanicalForce(Vec) and mechErrorNorm
   */
  void addNonconservativeForcing(double timeStep);

  /**
   * @brief Compute geometric forces, including
   * - spontaneous curvature force
   * - deviatoric curvature force
   * - area difference force
   * - capillary (surface tension) force
   * - osmotic force
   * - line capillary (line tension) force
   * - adsorption (area expansion) force
   * - aggregation (area expansion) force
   * - entropy (area expansion) force
   */
  void computeGeometricForces();
  void computeGeometricForces(size_t i);
  void computeGeometricForces(gcs::Vertex &v);

  /**
   * @brief Compute regularization pressure component of the system
   */
  void computeSpringForces();

  /**
   * @brief Compute Self Avoidance force
   */
  void computeSelfAvoidanceForce();

  /**
   * @brief Compute external force component of the system
   */
  EigenVectorX3dr prescribeExternalForce();

  /**
   * @brief Compute DPD forces of the system
   */
  void computeDPDForces(double dt);

  /**
   * @brief Compute damping forces of the system
   */
  gc::VertexData<gc::Vector3> computeDampingForce();

  /**
   * @brief Compute chemical potential of the system, including
   * - spontaneous curvature potential
   * - adsorption potential
   * - aggregation potential
   * - entropy potential
   * - deviatoric curvature potential
   * - dirichlet potential
   * - interior penalty potential
   */
  void computeChemicalPotentials();

  void computeNChemicalPotentials(int j);

  /**
   * @brief Compute in plane flux form on edge
   */
  EigenVectorX1d computeInPlaneFluxForm(EigenVectorX1d &chemicalPotential);

  EigenVectorX1d computeInPlaneFluxForm(gcs::VertexData<double> &pDensity, EigenVectorX1d &chemicalPotential);

  // ==========================================================
  // ================        energy.cpp      ==================
  // ==========================================================
  /**
   * @brief Compute spontaneous curvature energy
   */
  void computeSpontaneousCurvatureEnergy();

  /**
   * @brief Compute deviatoric curvature energy
   */
  void computeGaussianCurvatureEnergy();

  /**
   * @brief Compute deviatoric curvature energy
   */
  void computeDeviatoricCurvatureEnergy();

  /**
   * @brief Compute area difference energy
   */
  void computeAreaDifferenceEnergy();

  /**
   * @brief Compute surface energy
   */
  void computeSurfaceEnergy();

  /**
   * @brief Compute pressure work
   */
  void computePressureEnergy();

  /**
   * @brief Compute adsorption energy
   */
  void computeAdsorptionEnergy();

  /**
   * @brief Compute aggregation energy
   */
  void computeAggregationEnergy();

  /**
   * @brief Compute entropy penalty
   */
  void computeEntropyEnergy();

  /**
   * @brief Compute protein interior penalty
   */
  void computeProteinInteriorPenalty();

  /**
   * @brief Compute Dirichlet energy
   */
  void computeDirichletEnergy();

  /**
   * @brief Compute self-avoidance energy
   */
  void computeSelfAvoidanceEnergy();

  /**
   * @brief Compute edge spring energy
   */
  void computeEdgeSpringEnergy();

  /**
   * @brief Compute face spring energy
   */
  void computeFaceSpringEnergy();

  /**
   * @brief Compute LCR spring energy
   */
  void computeLcrSpringEnergy();

  /**
   * @brief Compute external work
   */
  double computeExternalWork(double currentTime, double dt);

  /**
   * @brief Compute kinetic energy
   */
  double computeKineticEnergy();

  /**
   * @brief Compute potential energy
   */
  double computePotentialEnergy();

  /**
   * @brief compute total energy
   */
  double computeTotalEnergy();

  /**
   * @brief Intermediate function to integrate the power
   */
  double computeIntegratedPower(double dt);
  double computeIntegratedPower(double dt, EigenVectorX3dr &&velocity);

  // ==========================================================
  // =============        regularization.cpp    ===============
  // ==========================================================
  /**
   * @brief Perform mesh mutation
   *
   * @param nMutation Iterations of mutations
   */
  void mutateMesh(size_t nMutation = 1);

  /**
   * @brief Apply vertex shift by moving the vertices chosen for integration to
   * the Barycenter of the it neighbors
   */
  void vertexShift();

  /**
   * @brief One pass edge flipping if not Delaunay
   */
  bool edgeFlip();
  /**
   * @brief Flip edges using a queue if not Delaunay
   *
   */
  void edgeFlipQueued();

  inline gcs::Vertex splitEdge(gcs::Edge e) {
    gcs::Vertex vertex1 = e.firstVertex(), vertex2 = e.secondVertex();
    gc::Vector3 vertex1Pos = geometry.vpg->vertexPositions[vertex1];
    gc::Vector3 vertex2Pos = geometry.vpg->vertexPositions[vertex2];
    gc::Vector3 vertex1Vel = velocity[vertex1];
    gc::Vector3 vertex2Vel = velocity[vertex2];
    double vertex1GeoDist = geometry.geodesicDistance[vertex1];
    double vertex2GeoDist = geometry.geodesicDistance[vertex2];
    double vertex1Phi = proteinDensity[vertex1];
    double vertex2Phi = proteinDensity[vertex2];
    double vertex1Area = geometry.vpg->vertexDualAreas[vertex1];
    double vertex2Area = geometry.vpg->vertexDualAreas[vertex2];
    std::vector<double> vertex1Phis;
    std::vector<double> vertex2Phis;
    for (int j = 0; j < pDensities.size(); j++){
      vertex1Phis.push_back(pDensities[j][vertex1]);
      vertex2Phis.push_back(pDensities[j][vertex2]);
    }

    // bool vertex1PointTracker = geometry.notableVertex[vertex1];
    // bool vertex2PointTracker = geometry.notableVertex[vertex2];

    // split the edge
    gcs::Vertex newVertex = geometry.mesh->splitEdgeTriangular(e).vertex();

    // update quantities
    // Note: think about conservation of energy, momentum and angular
    // momentum
    // averageData(geometry.vpg->inputVertexPositions, vertex1, vertex2,
    // newVertex); averageData(velocity, vertex1, vertex2, newVertex);
    // averageData(geometry.geodesicDistance, vertex1, vertex2, newVertex);
    // averageData(proteinDensity, vertex1, vertex2, newVertex);
    geometry.vpg->vertexPositions[newVertex] = 0.5 * (vertex1Pos + vertex2Pos);
    velocity[newVertex] = 0.5 * (vertex1Vel + vertex2Vel);
    geometry.geodesicDistance[newVertex] =
        0.5 * (vertex1GeoDist + vertex2GeoDist);
    proteinDensity[newVertex] = 0.5 * (vertex1Phi + vertex2Phi);
    for (int j = 0; j < pDensities.size(); j++){
      if (pParameters[j].conserve)
        pDensities[j][newVertex] =
            (vertex1Phis[j] *
                (vertex1Area - geometry.vpg->vertexDualAreas[vertex1]) +
            vertex2Phis[j] *
                (vertex2Area - geometry.vpg->vertexDualAreas[vertex2])) /
            geometry.vpg->vertexDualAreas[newVertex];
      else
        pDensities[j][newVertex] = 0.5 * (vertex1Phis[j] + vertex2Phis[j]);
    }
    geometry.notableVertex[newVertex] = false;
    forces.forceMask[newVertex] = gc::Vector3{1, 1, 1};
    return newVertex;
  }

  /**
   * @brief Perform edge collapse
   *
   * @warning if e is a part of a pinch triangle than the new vertex can be
   * null.
   *
   * @param e Edge to collapse
   * @return gcs::Vertex Resultant vertex. Warning: if the edge collapse would
   * cause issues the result can be null
   */
  inline gcs::Vertex collapseEdge(gcs::Edge e) {
    // First collect data
    gcs::Vertex vertex1 = e.firstVertex(), vertex2 = e.secondVertex();
    gc::Vector3 vertex1ForceMask = forces.forceMask[vertex1];
    gc::Vector3 vertex2ForceMask = forces.forceMask[vertex2];
    gc::Vector3 vertex1Pos = geometry.vpg->vertexPositions[vertex1];
    gc::Vector3 vertex2Pos = geometry.vpg->vertexPositions[vertex2];
    gc::Vector3 vertex1Vel = velocity[vertex1];
    gc::Vector3 vertex2Vel = velocity[vertex2];
    double vertex1GeoDist = geometry.geodesicDistance[vertex1];
    double vertex2GeoDist = geometry.geodesicDistance[vertex2];
    double vertex1Phi = proteinDensity[vertex1];
    double vertex2Phi = proteinDensity[vertex2];
    std::vector<double> vertex1Phis;
    std::vector<double> vertex2Phis;
    for (int j = 0; j < pDensities.size(); j++) {
      vertex1Phis.push_back(pDensities[j][vertex1]);
      vertex2Phis.push_back(pDensities[j][vertex2]);
    }
    std::vector<double> totalProteins;
    for (int j = 0; j < pDensities.size(); j++) {
      totalProteins.push_back(
          pDensities[j][vertex1] * geometry.vpg->vertexDualAreas[vertex1] +
          pDensities[j][vertex2] * geometry.vpg->vertexDualAreas[vertex2]);
    }
    bool vertex1PointTracker = geometry.notableVertex[vertex1];
    bool vertex2PointTracker = geometry.notableVertex[vertex2];

    // collapse the edge
    gcs::Vertex newVertex = geometry.mesh->collapseEdgeTriangular(e);

    // vertex maybe void if edge was on a pinch triangle
    if (newVertex.getIndex() != gc::INVALID_IND) {

      // update quantities
      // Note: think about conservation of energy, momentum and angular
      // momentum
      geometry.vpg->vertexPositions[newVertex] =
          ((gc::sum(vertex1ForceMask) < 2.5) || vertex1PointTracker)
              ? vertex1Pos
          : ((gc::sum(vertex2ForceMask) < 2.5) || vertex2PointTracker)
              ? vertex2Pos
              : (vertex1Pos + vertex2Pos) / 2;
      // averageData(velocity, vertex1, vertex2, newVertex);
      // averageData(geometry.geodesicDistance, vertex1, vertex2, newVertex);
      // averageData(proteinDensity, vertex1, vertex2, newVertex);
      velocity[newVertex] = 0.5 * (vertex1Vel + vertex2Vel);
      geometry.geodesicDistance[newVertex] =
          0.5 * (vertex1GeoDist + vertex2GeoDist);
      proteinDensity[newVertex] = 0.5 * (vertex1Phi + vertex2Phi);
      for (int j = 0; j < pDensities.size(); j++){
        if (pParameters[j].conserve)
          pDensities[j][newVertex] =
              totalProteins[j] / geometry.vpg->vertexDualAreas[newVertex];
        else
          pDensities[j][newVertex] = 0.5 * (vertex1Phis[j] + vertex2Phis[j]);
      }
      geometry.notableVertex[newVertex] =
          vertex1PointTracker || vertex2PointTracker;
    }
    return newVertex;
  }

  bool ifFoldover(const gcs::Edge e, const double angle = 0.5) const;

  /**
   * @brief Split/collapse mesh mutations without cache
   *
   * @return true
   * @return false
   */
  bool processSplitCollapse();

  /**
   * @brief Process split/collapse mutations with a queue
   *
   * @return true
   * @return false
   */
  bool processSplitCollapseQueued();

  /**
   * @brief global smoothing after mutation of the mesh
   * @param initStep initial guess of time step
   * @param target target reduce of force norm
   * @param maxIteration maximum number of iteration
   */
  Eigen::Matrix<bool, Eigen::Dynamic, 1>
  smoothenMesh(double initStep = 0.01, double target = 0.7,
               size_t maxIteration = 1000);
  /**
   * @brief pointwise smoothing after mutation of the mesh
   */
  void localSmoothing(const gcs::Vertex &v, std::size_t num = 10,
                      double stepSize = 0.01);
  void localSmoothing(const gcs::Halfedge &he, std::size_t num = 10,
                      double stepSize = 0.01);

  /**
   * @brief global update of quantities after mutation of the mesh
   */
  void globalUpdateAfterMutation();

  // ==========================================================
  // =============          misc.cpp            ===============
  // ==========================================================
  /**
   * @brief test conservative force computation by validating energy decrease
   * @return
   */
  bool testConservativeForcing(const double timeStep);

  /**
   * @brief backtrace energy increase from the system perturbation
   * @return
   */
  void backtraceEnergyGrowth(const double timeStep,
                             const Energy previousEnergy);

  /**
   * @brief Check finiteness of forcing and energy
   * @return whether the system is finite
   */
  bool checkFiniteness();

  /**
   * @brief testing of random number generator pcg
   *
   */
  // void check_pcg();

  /**
   * @brief prescribe mask based on geodesic disk
   */
  void prescribeGeodesicMasks();
};
} // namespace solver
} // namespace mem3dg
