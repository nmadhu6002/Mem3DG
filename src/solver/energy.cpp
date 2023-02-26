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

#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::computeSpontaneousCurvatureEnergy() {

  // no particular reason, just experiment
  // E.BE = 0;
  // for (gcs::Vertex v : mesh->vertices()) {
  //   E.BE += Kb[v] *
  //           (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
  //           (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
  //           vpg->vertexDualArea(v);
  // }

  // Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H.raw() - H0.raw();
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference;

  Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference =
      abs(vpg->vertexMeanCurvatures.raw().array() /
              vpg->vertexDualAreas.raw().array() -
          H0.raw().array());
  energy.spontaneousCurvatureEnergy =
      (Kb.raw().array() * vpg->vertexDualAreas.raw().array() *
       H_difference.array().square())
          .sum();

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeDeviatoricCurvatureEnergy() {
  energy.deviatoricCurvatureEnergy =
      (Kd.raw().array() * vpg->vertexGaussianCurvatures.raw().array().square())
          .sum();
  // (Kd.raw().array() * (vpg->vertexMeanCurvatures.raw().array().square() /
  //                          vpg->vertexDualAreas.raw().array() -
  //                      vpg->vertexGaussianCurvatures.raw().array()))
  //     .sum();
}

void System::computeAreaDifferenceEnergy() {
  double K =
      0.25 * parameters.bending.alpha * parameters.bending.Kb * constants::PI;
  energy.areaDifferenceEnergy =
      K / parameters.bending.D / parameters.bending.D / surfaceArea *
      pow(2 * parameters.bending.D * vpg->vertexMeanCurvatures.raw().sum() -
              parameters.bending.dA0,
          2);
}

void System::computeSurfaceEnergy() {
  // cotan laplacian normal is exact for area variation
  double A_difference = surfaceArea - parameters.tension.At;
  energy.surfaceEnergy =
      parameters.tension.isConstantSurfaceTension
          ? forces.surfaceTension * surfaceArea
          : forces.surfaceTension * A_difference / 2 +
                parameters.tension.lambdaSG * A_difference / 2;
}

void System::computeEdgeSpringEnergy() {
  energy.edgeSpringEnergy =
      0.5 * parameters.spring.Kse *
      ((vpg->edgeLengths.raw() - refVpg->edgeLengths.raw()).array() /
       refVpg->edgeLengths.raw().array())
          .square()
          .sum();
}

void System::computeFaceSpringEnergy() {
  energy.faceSpringEnergy =
      0.5 * parameters.spring.Ksl *
      ((vpg->faceAreas.raw() - refVpg->faceAreas.raw()).array() /
       refVpg->faceAreas.raw().array())
          .square()
          .sum();
}

void System::computeLcrSpringEnergy() {
  double E = 0;
  for (std::size_t i = 0; i < mesh->nEdges(); ++i) {
    gc::Edge e{mesh->edge(i)};
    gc::Halfedge he = e.halfedge();
    double lcr = computeLengthCrossRatio(*vpg, he);
    E += pow((lcr - refLcrs[he]) / refLcrs[he], 2);
  }
  E *= 0.5 * parameters.spring.Kst;
  energy.lcrSpringEnergy = E;
}

void System::computePressureEnergy() {
  // Note: area weighted normal is exact volume variation
  if (parameters.osmotic.isPreferredVolume) {
    double V_difference = volume - parameters.osmotic.Vt;
    energy.pressureEnergy = -forces.osmoticPressure * V_difference / 2 +
                            parameters.osmotic.lambdaV * V_difference / 2;
  } else if (parameters.osmotic.isConstantOsmoticPressure) {
    energy.pressureEnergy = -forces.osmoticPressure * volume;
  } else {
    double ratio = parameters.osmotic.cam * volume / parameters.osmotic.n;
    energy.pressureEnergy = mem3dg::constants::i * mem3dg::constants::R *
                            parameters.temperature * parameters.osmotic.n *
                            (ratio - log(ratio) - 1);
  }
}

// void System::computeAdsorptionEnergy() {
//   energy.adsorptionEnergy =
//       parameters.adsorption.epsilon * (proteinDensity.raw().array()).sum();
// }

// void System::computeAggregationEnergy() {
//   energy.aggregationEnergy =
//       parameters.aggregation.chi *
//       (proteinDensity.raw().array() * proteinDensity.raw().array()).sum();
// }

void System::computeAdsorptionEnergy(bool protein, bool protein2) {
  // std::cout << "computeAdsorptionEnergy" << std::endl;
  if (protein)
    energy.adsorptionEnergy =
        parameters.adsorption.epsilon *
        (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array()).sum();
  if (protein2)
    energy.adsorption2Energy =
        parameters.adsorption.epsilon2 *
        (vpg->vertexDualAreas.raw().array() * protein2Density.raw().array()).sum();
}

void System::computeAggregationEnergy(bool protein, bool protein2) {
  // std::cout << "computeAggregationEnergy" << std::endl;
  if (protein)
    energy.aggregationEnergy =
        parameters.aggregation.chi *
        (vpg->vertexDualAreas.raw().array() *
        ((2 * proteinDensity.raw().array() - 1).square() - 1).square())
            .sum();
  if (protein2)
    energy.aggregation2Energy =
        parameters.aggregation.chi2 *
        (vpg->vertexDualAreas.raw().array() *
        ((2 * protein2Density.raw().array() - 1).square() - 1).square())
            .sum();
  // energy.aggregationEnergy =
  //     parameters.aggregation.chi *
  //     (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array() *
  //      proteinDensity.raw().array())
  //         .sum();
}

void System::computeEntropyEnergy(bool protein, bool protein2) {
  // std::cout << "computeEntropyEnergy" << std::endl;
  if (protein)
    energy.entropyEnergy =
        parameters.entropy.xi *
        (vpg->vertexDualAreas.raw().array() *
        (proteinDensity.raw().array().log() * proteinDensity.raw().array() +
          (1 - proteinDensity.raw().array()).log() *
              (1 - proteinDensity.raw().array())))
            .sum();
  if (protein2)
    energy.entropy2Energy =
        parameters.entropy.xi2 *
        (vpg->vertexDualAreas.raw().array() *
        (protein2Density.raw().array().log() * protein2Density.raw().array() +
          (1 - protein2Density.raw().array()).log() *
              (1 - protein2Density.raw().array())))
            .sum();
  // energy.entropyEnergy =
  //     parameters.entropy.xi *
  //     (vpg->vertexDualAreas.raw().array() *
  //      (proteinDensity.raw().array().log() - 1) *
  //      proteinDensity.raw().array())
  //         .sum();
}

void System::computeProteinInteriorPenalty(bool protein, bool protein2) {
  // std::cout << "computeProteinInteriorPenalty" << std::endl;
  // interior method to constrain protein density to remain from 0 to 1
  if (protein)
    energy.proteinInteriorPenalty =
        -parameters.protein.proteinInteriorPenalty *
        ((proteinDensity.raw().array()).log().sum() +
        (1 - proteinDensity.raw().array()).log().sum());
  if (protein2)
    energy.protein2InteriorPenalty =
        -parameters.protein2.proteinInteriorPenalty *
        ((protein2Density.raw().array()).log().sum() +
        (1 - protein2Density.raw().array()).log().sum());
}

void System::computeSelfAvoidanceEnergy() {
  const double d0 = parameters.selfAvoidance.d;
  const double mu = parameters.selfAvoidance.mu;
  const double n = parameters.selfAvoidance.n;
  double e = 0.0;
  projectedCollideTime = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex vi{mesh->vertex(i)};
    gc::VertexData<bool> neighborList(*mesh, false);
    meshProcessor.meshMutator.markVertices(neighborList, vi, n);
    for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
      if (neighborList[j])
        continue;
      gc::Vertex vj{mesh->vertex(j)};

      // double penalty = mu * vpg->vertexDualAreas[vi] * proteinDensity[vi] *
      //                  vpg->vertexDualAreas[vj] * proteinDensity[vj];
      double penalty = mu * proteinDensity[vi] * proteinDensity[vj];
      // double penalty = mu;
      // double penalty = mu * vpg->vertexDualAreas[vi] *
      // vpg->vertexDualAreas[vj];

      gc::Vector3 r =
          vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
      double distance = gc::norm(r) - d0;
      double collideTime = distance / gc::dot(velocity[vi] - velocity[vj], r);
      if (collideTime < projectedCollideTime &&
          gc::dot(velocity[vi] - velocity[vj], r) > 0)
        projectedCollideTime = collideTime;
      // e -= penalty * log(distance);
      e += penalty / distance;
    }
  }
  if (projectedCollideTime == std::numeric_limits<double>::max())
    projectedCollideTime = 0;
  energy.selfAvoidancePenalty = e;
}

void System::computeDirichletEnergy(bool protein, bool protein2) {
  // std::cout << "computeDirichletEnergy" << std::endl;
  if (false) {
    mem3dg_runtime_error("computeDirichletEnergy: out of date implementation, "
                         "shouldn't be called!");
    // scale the dH0 such that it is integrated over the edge
    // this is under the case where the resolution is low, WIP
    // auto dH0 = vpg->edgeLengths.raw().array() *  ((vpg->d0 *
    // H0.raw()).cwiseAbs()).array(); auto dH0 = (vpg->d0 *
    // H0.raw()).cwiseAbs();
    // E.dE = (vpg->hodge1Inverse * F.lineTension.raw()).sum();
  }

  // explicit dirichlet energy
  if (protein){
    energy.dirichletEnergy = 0;
    for (gcs::Face f : mesh->faces()) {
      energy.dirichletEnergy += 0.5 * parameters.dirichlet.eta *
                                proteinDensityGradient[f].norm2() *
                                vpg->faceAreas[f];
    }
  }
  if (protein2){
    energy.dirichlet2Energy = 0;
    for (gcs::Face f : mesh->faces()) {
      energy.dirichlet2Energy += 0.5 * parameters.dirichlet.eta2 *
                                 protein2DensityGradient[f].norm2() *
                                 vpg->faceAreas[f];
    }
  }

  // alternative dirichlet energy after integration by part
  // E.dE = 0.5 * P.eta * proteinDensity.raw().transpose() *
  // vpg->cotanLaplacian
  // *
  //        proteinDensity.raw();
}

double System::computePotentialEnergy() {
  // std::cout << "computePotentialEnergy" << std::endl;
  // fundamental internal potential energy
  energy.spontaneousCurvatureEnergy = 0;
  energy.deviatoricCurvatureEnergy = 0;
  energy.areaDifferenceEnergy = 0;
  energy.surfaceEnergy = 0;
  energy.pressureEnergy = 0;
  energy.adsorptionEnergy = 0;
  energy.aggregationEnergy = 0;
  energy.entropyEnergy = 0;
  energy.dirichletEnergy = 0;
  energy.adsorption2Energy = 0;
  energy.aggregation2Energy = 0;
  energy.entropy2Energy = 0;
  energy.dirichlet2Energy = 0;
  energy.selfAvoidancePenalty = 0;
  energy.proteinInteriorPenalty = 0;
  energy.protein2InteriorPenalty = 0;
  energy.edgeSpringEnergy = 0;
  energy.faceSpringEnergy = 0;
  energy.lcrSpringEnergy = 0;

  computeSpontaneousCurvatureEnergy();
  computeDeviatoricCurvatureEnergy();

  // optional internal potential energy
  if (parameters.bending.alpha != 0)
    computeAreaDifferenceEnergy();
  if (parameters.tension.Ksg != 0)
    computeSurfaceEnergy();
  if (parameters.osmotic.Kv != 0)
    computePressureEnergy();
  computeAdsorptionEnergy(parameters.adsorption.epsilon != 0,
                          parameters.adsorption.epsilon2 != 0);
  computeAggregationEnergy(parameters.aggregation.chi != 0,
                           parameters.aggregation.chi2 != 0);
  computeEntropyEnergy(parameters.entropy.xi != 0,
                       parameters.entropy.xi2 != 0);
  computeDirichletEnergy(parameters.dirichlet.eta != 0,
                         parameters.dirichlet.eta2 != 0);
  if (parameters.selfAvoidance.mu != 0) {
    computeSelfAvoidanceEnergy();
  }
  if (parameters.spring.Kse != 0) {
    computeEdgeSpringEnergy();
  }
  if (parameters.spring.Ksl != 0) {
    computeFaceSpringEnergy();
  }
  if (parameters.spring.Kst != 0) {
    computeLcrSpringEnergy();
  }
  computeProteinInteriorPenalty(
      parameters.variation.isProteinVariation &&
          parameters.protein.proteinInteriorPenalty != 0,
      parameters.variation.isProtein2Variation &&
          parameters.protein2.proteinInteriorPenalty != 0);

  // summerize internal potential energy
  energy.potentialEnergy =
      energy.spontaneousCurvatureEnergy + energy.deviatoricCurvatureEnergy +
      energy.areaDifferenceEnergy + energy.surfaceEnergy +
      energy.pressureEnergy + energy.adsorptionEnergy + energy.dirichletEnergy +
      energy.aggregationEnergy + energy.entropyEnergy + energy.adsorption2Energy + 
      energy.dirichlet2Energy + energy.aggregation2Energy + energy.entropy2Energy +
      energy.selfAvoidancePenalty + energy.proteinInteriorPenalty + 
      energy.proteinInteriorPenalty + energy.edgeSpringEnergy + energy.faceSpringEnergy +
      energy.lcrSpringEnergy;
  return energy.potentialEnergy;
}

double System::computeIntegratedPower(double dt) {
  prescribeExternalForce();
  return dt * rowwiseDotProduct(toMatrix(forces.externalForceVec),
                                toMatrix(velocity))
                  .sum();
}
double System::computeIntegratedPower(double dt, EigenVectorX3dr &&velocity) {
  prescribeExternalForce();
  return dt *
         rowwiseDotProduct(toMatrix(forces.externalForceVec), velocity).sum();
}

double System::computeExternalWork(double currentTime, double dt) {
  energy.time = currentTime;
  energy.externalWork += computeIntegratedPower(dt);
  return energy.externalWork;
}

double System::computeKineticEnergy() {
  energy.kineticEnergy =
      0.5 * Eigen::square(toMatrix(velocity).array()).matrix().sum();
  // energy.kineticEnergy =
  //     0.5 * rowwiseScalarProduct(toMatrix(vpg->vertexDualAreas).array(),
  //                                Eigen::square(toMatrix(velocity).array()))
  //               .sum();
  return energy.kineticEnergy;
}

double System::computeTotalEnergy() {
  // std::cout << "computeTotalEnergy" << std::endl;
  computePotentialEnergy();
  computeKineticEnergy();
  if (time == energy.time) {
    energy.totalEnergy =
        energy.kineticEnergy + energy.potentialEnergy - energy.externalWork;
  } else if (!parameters.external.isActivated) {
    energy.totalEnergy = energy.kineticEnergy + energy.potentialEnergy;
  } else {
    mem3dg_runtime_error("energy.externalWork not updated!")
  }
  return energy.totalEnergy;
}

void System::computeFaceTangentialDerivative(
    gcs::VertexData<double> &quantities, gcs::FaceData<gc::Vector3> &gradient) {
  if ((quantities.raw().array() == quantities.raw()[0]).all()) {
    gradient.fill({0, 0, 0});
  } else {
    for (gcs::Face f : mesh->faces()) {
      gc::Vector3 normal = vpg->faceNormals[f];
      gc::Vector3 gradientVec{0, 0, 0};
      for (gcs::Halfedge he : f.adjacentHalfedges()) {
        gradientVec += quantities[he.next().tipVertex()] *
                       gc::cross(normal, vecFromHalfedge(he, *vpg));
      }
      gradient[f] = gradientVec / 2 / vpg->faceAreas[f];
    }
  }
}

double System::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                       gcs::Halfedge &he) const {
  if (he.edge().isBoundary()) {
    return 1;
  } else {
    gcs::Edge lj = he.next().edge();
    gcs::Edge ki = he.twin().next().edge();
    gcs::Edge il = he.next().next().edge();
    gcs::Edge jk = he.twin().next().next().edge();
    return vpg.edgeLengths[il] * vpg.edgeLengths[jk] / vpg.edgeLengths[ki] /
           vpg.edgeLengths[lj];
  }
}

} // namespace solver
} // namespace mem3dg
