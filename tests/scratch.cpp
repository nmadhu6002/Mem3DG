#include <iostream>

#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>

// We are writing 2D data, a 6 x 12 grid
constexpr int nx = 6;
constexpr int ny = 12;

// Return this in event of a problem
constexpr int nc_err = 2;

int main() {
  namespace gc = ::geometrycentral;
  namespace gcs = ::geometrycentral::surface;

  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  ddgsolver::tetrahedron(coords, polygons);

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();

  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);

  auto file = ddgsolver::trajfile::TrajFile::newFile("test.nc", *ptrMesh, true);

  // file.writeTime(file.getNextFrameIndex(), 1);
  // file.writeTime(file.getNextFrameIndex(), 2);

  // file.writeCoords(
  //     0, ddgsolver::EigenMap<double, 3>(ptrVpg->inputVertexPositions));
  // file.writeCoords(
  //     3, ddgsolver::EigenMap<double, 3>(ptrVpg->inputVertexPositions));

  // double x, y;
  // ddgsolver::TrajFile::EigenVector vec1, vec2;

  // std::tie(x, vec1) = file.getTimeAndCoords(0);
  // std::cout << "Time " << x << std::endl << vec1 << std::endl;
  
  // auto file2 = ddgsolver::TrajFile::openReadOnly("test.nc");
  // std::tie(y, vec2) = file2.getTimeAndCoords(1);
  // std::cout << "Time " << y << std::endl << vec2 << std::endl;

  std::cout << "EOF" << std::endl;
  return 0;
}
