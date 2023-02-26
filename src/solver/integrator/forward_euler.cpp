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

#include <Eigen/Core>
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool Euler::integrate() {
  std::cout << "integrate" << std::endl;
  if (ifDisableIntegrate)
    mem3dg_runtime_error("integrate() is disabled for current construction!");

  signal(SIGINT, signalHandler);

  double initialTime = system.time, lastUpdateGeodesics = system.time,
         lastProcessMesh = system.time, lastComputeAvoidingForce = system.time,
         lastSave = system.time;

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (ifOutputTrajFile) {
    createMutableNetcdfFile(isContinuation);
    if (ifPrintToConsole)
      std::cout << "Initialized NetCDF file at "
                << outputDirectory + "/" + trajFileName << std::endl;
  }
#endif

  // time integration loop
  const double avoidStrength = system.parameters.selfAvoidance.mu;
  unsigned count = 0;
  for (;;) {
    std::cout << count++  << std::endl;
    // turn on/off self-avoidance; outside status-march-cycle; before savedata
    // to write selfAvoidance
    if (avoidStrength != 0) {
      if ((system.time - lastComputeAvoidingForce) >
              system.parameters.selfAvoidance.p * system.projectedCollideTime ||
          system.time - lastSave >= savePeriod || system.time == initialTime ||
          EXIT) {
        lastComputeAvoidingForce = system.time;
        system.parameters.selfAvoidance.mu = avoidStrength;
        if (ifPrintToConsole) {
          std::cout << "computing avoiding force at "
                    << "t = " << system.time << std::endl;
          std::cout << "projected collision is " << system.projectedCollideTime
                    << std::endl;
          std::cout << "time step is " << timeStep << std::endl;
        }
      } else {
        system.parameters.selfAvoidance.mu = 0;
      }
    }

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info; save data before exit
    if (system.time - lastSave >= savePeriod || system.time == initialTime ||
        EXIT) {
      lastSave = system.time;
      saveData(ifOutputTrajFile, ifOutputMeshFile, ifPrintToConsole);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // Process mesh every tProcessMesh period
    if (system.time - lastProcessMesh > (processMeshPeriod * timeStep)) {
      lastProcessMesh = system.time;
      system.mutateMesh();
      system.updateConfigurations();
      system.refVpg = system.vpg->copy();
      system.updateReferenceConfigurations();
      if (system.parameters.point.isFloatVertex)
        system.findFloatCenter(
            3 * system.vpg->edgeLength(
                    system.center.nearestVertex().halfedge().edge()));
    }

    // update geodesics every tUpdateGeodesics period
    if (system.time - lastUpdateGeodesics >
        (updateGeodesicsPeriod * timeStep)) {
      lastUpdateGeodesics = system.time;
      if (system.parameters.point.isFloatVertex)
        system.findFloatCenter(
            3 * system.vpg->edgeLength(
                    system.center.nearestVertex().halfedge().edge()));
      system.geodesicDistance.raw() = system.computeGeodesicDistance();
      if (system.parameters.protein.ifPrescribe)
        system.prescribeGeodesicProteinDensityDistribution();
      system.updateConfigurations();
    }

    // step forward
    if (system.time == lastProcessMesh || system.time == lastUpdateGeodesics) {
      system.time += 1e-10 * characteristicTimeStep;
    } else {
      march();
    }
  }

#ifdef MEM3DG_WITH_NETCDF
  if (ifOutputTrajFile) {
    closeMutableNetcdfFile();
    if (ifPrintToConsole)
      std::cout << "Closed NetCDF file" << std::endl;
  }
#endif

  // return if optimization is sucessful
  if (!SUCCESS && ifOutputTrajFile) {
    std::string filePath = outputDirectory;
    filePath.append("/");
    filePath.append(trajFileName);
  }

  return SUCCESS;
}

void Euler::checkParameters() {
  if (system.parameters.dpd.gamma != 0) {
    mem3dg_runtime_error("DPD has to be turned off for euler integration!");
  }
  if (system.parameters.damping != 0) {
    mem3dg_runtime_error("Damping to be 0 for euler integration!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      mem3dg_runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
}

void Euler::status() {
  std::cout << "status" << std::endl;
  // compute summerized forces
  system.computeConservativeForcing();
  system.addNonconservativeForcing(timeStep);

  // exit if under error tolerance
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance && system.chem2ErrorNorm < tolerance) {
    if (ifPrintToConsole)
      std::cout << "\nError norm smaller than tolerance." << std::endl;
    EXIT = true;
  }

  // exit if reached time
  if (system.time > totalTime) {
    if (ifPrintToConsole)
      std::cout << "\nReached time." << std::endl;
    EXIT = true;
  }

  // compute the free energy of the system
  if (system.parameters.external.isActivated)
    system.computeExternalWork(system.time, timeStep);
  system.computeTotalEnergy();

  // check finiteness
  if (!std::isfinite(timeStep) || !system.checkFiniteness()) {
    EXIT = true;
    SUCCESS = false;
    if (!std::isfinite(timeStep))
      mem3dg_runtime_message("time step is not finite!");
  }
}

void Euler::march() {
  std::cout << "march" << std::endl;
  // compute velocity, which are independent of time
  system.velocity = system.forces.mechanicalForceVec;
  system.mechErrorNorm = (toMatrix(system.velocity).array() *
                          toMatrix(system.forces.mechanicalForceVec).array())
                             .sum();
  if (system.parameters.variation.isProteinVariation) {
    if (system.parameters.variation.isProteinConservation) {
      system.proteinRateOfChange.raw() =
          system.parameters.proteinMobility * system.vpg->hodge0Inverse *
          system.vpg->d0.transpose() *
          system.computeInPlaneFluxForm(system.forces.chemicalPotential.raw(), system.proteinDensity);
    } else {
      system.proteinRateOfChange = system.parameters.proteinMobility *
                                   system.forces.chemicalPotential /
                                   system.vpg->vertexDualAreas;
    }
    system.chemErrorNorm = (system.proteinRateOfChange.raw().array() *
                            system.forces.chemicalPotential.raw().array())
                               .sum();
  }

  if (system.parameters.variation.isProtein2Variation) {
    if (system.parameters.variation.isProtein2Conservation) {
      system.protein2RateOfChange.raw() =
          system.parameters.protein2Mobility * system.vpg->hodge0Inverse *
          system.vpg->d0.transpose() *
          system.computeInPlaneFluxForm(system.forces.chemical2Potential.raw(), system.protein2Density);
    } else {
      system.protein2RateOfChange = system.parameters.protein2Mobility *
                                   system.forces.chemical2Potential /
                                   system.vpg->vertexDualAreas;
    }
    system.chem2ErrorNorm = (system.protein2RateOfChange.raw().array() *
                            system.forces.chemical2Potential.raw().array())
                               .sum();
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (ifAdaptiveStep) {
    characteristicTimeStep = getAdaptiveCharacteristicTimeStep();
  }

  // backtracking to obtain stable time step
  if (isBacktrack) {
    double timeStep_mech = std::numeric_limits<double>::max(),
           timeStep_chem = std::numeric_limits<double>::max(),
           timeStep_chem2 = std::numeric_limits<double>::max();
    if (system.parameters.variation.isShapeVariation)
      timeStep_mech = mechanicalBacktrack(toMatrix(system.velocity), rho, c1);
    if (system.parameters.variation.isProteinVariation)
      timeStep_chem =
          chemicalBacktrack(system.proteinRateOfChange.raw(), rho, c1);
    if (system.parameters.variation.isProtein2Variation)
      timeStep_chem2 =
          chemical2Backtrack(system.protein2RateOfChange.raw(), rho, c1);
    double temp_timeStep =(timeStep_chem < timeStep_mech) ? timeStep_chem : timeStep_mech;
    timeStep = (timeStep_chem2 < temp_timeStep) ? timeStep_chem2 : temp_timeStep;
  } else {
    timeStep = characteristicTimeStep;
  }
  system.vpg->inputVertexPositions += system.velocity * timeStep;
  system.proteinDensity += system.proteinRateOfChange * timeStep;
  system.protein2Density += system.protein2RateOfChange * timeStep;
  system.time += timeStep;

  // recompute cached values
  system.updateConfigurations();
  std::cout << "march end" << std::endl;
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
