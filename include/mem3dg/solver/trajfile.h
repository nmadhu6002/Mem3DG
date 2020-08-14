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

/**
 * @file  trajfile.h
 * @brief Netcdf trajectory output support
 *
 */

#pragma once

#ifdef MEM3DG_WITH_NETCDF

#include <exception>
#include <iostream>
#include <sstream>
#include <vector>

#include <Eigen/Core>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/typetraits.h"

namespace ddgsolver {

namespace trajfile {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// namespace nc = ::netCDF;
extern "C" {
#include <netcdf.h>
}

// DIMENSIONS NAMES
static const std::string POLYGON_ORDER_NAME = "polygon_dims";
/// Number of vertices per polygon
static const std::size_t POLYGON_ORDER = 3;

static const std::string NPOLYGONS_NAME = "npolygons";

static const std::string SPATIAL_DIMS_NAME = "spatial";
/// Number of spatial dimensions
static const std::size_t SPATIAL_DIMS = 3;
/// Name of frames
static const std::string FRAME_NAME = "frame";
/// nvertices
static const std::string NVERTICES_NAME = "nvertices";

/// Name of conventions
static const std::string CONVENTIONS_NAME = "Conventions";
/// Convention value
static const std::string CONVENTIONS_VALUE = "Mem3DG";
/// Conventions version
static const std::string CONVENTIONS_VERSION_NAME = "ConventionsVersion";
/// Conventions version value
static const std::string CONVENTIONS_VERSION_VALUE = "0.0.1";

/// Name of the units labels
static const std::string UNITS = "units";
/// Value for length units
static const std::string LENGTH_UNITS = "micrometers";
/// Value for time units
static const std::string TIME_UNITS = "picoseconds";

// Data/Variable block names
/// Name of time data
static const std::string TIME_DNAME = "time";
/// Name of coordinates data
static const std::string COORD_DNAME = "coordinates";
/// Name of the mesh topology data
static const std::string TOPOLOGY_DNAME = "topology";
/// Name of the velocity data
static const std::string VELOCITY_DNAME = "velocities";
/// Name of the mean curvature data
static const std::string MEAN_CURVATURE_DNAME = "mean_curvature";
/// Name of the gaussian curvature data
static const std::string GAUSSIAN_CURVATURE_DNAME = "gaussian_curvature";

namespace util {
// Handle errors by printing an error message and exiting with a non-zero
// status.
inline void ncCheck(int retCode, const char *file, int line) {
  if (retCode) {
    std::cout << retCode << std::endl;
    std::stringstream ss;
    ss << "In file " << file << " on line " << line << " "
       << nc_strerror(retCode);
    throw std::runtime_error(ss.str());
  }
}

template <int ndims>
inline void def_var(int ncid, const std::string &name, nc_type xtype,
                    const int (&dimidsp)[ndims], int *varidp, const char *file,
                    int line) {
  ncCheck(nc_def_var(ncid, name.c_str(), xtype, ndims, dimidsp, varidp), file,
          line);
}
} // namespace util

/**
 * @class TrajFile
 * @brief Trajectory interface to help with manipulating trajectories
 *
 */
class DLL_PUBLIC TrajFile {
private:
  TrajFile *operator=(const TrajFile &rhs) = delete;
  TrajFile(const TrajFile &rhs) = delete;

public:
  using EigenVector =
      Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>;

  TrajFile() = delete;

  TrajFile(TrajFile &&rhs) = default;

  /**
   * @brief Open a new file and populate it with the convention
   *
   * @param filename  Filename to save to
   * @param mesh      Mesh of interest
   * @param replace   Whether to replace an existing file or exit
   *
   * @return TrajFile helper object to manipulate the bound NetCDF file.
   */
  static TrajFile newFile(const std::string &filename, gcs::SurfaceMesh &mesh,
                          bool replace = false) {
    if (replace)
      return TrajFile(filename, mesh, NC_CLOBBER);
    else
      return TrajFile(filename, mesh, NC_NOCLOBBER);
  };

  /**
   * @brief Destructor frees the bound NcFile
   */
  ~TrajFile() {
    if (isOpen)
      util::ncCheck(nc_close(ncid), __FILE__, __LINE__);
  };

  inline bool isWriteable() { return writeable; };

private:
  /**
   * @brief Private constructor for opening an existing file.
   *
   * The metadata is checked for consistency with the convention. This
   * function should only be called with fMode set to `read` or `write`.
   *
   * @param filename Path to file of interest
   * @param fMode    Mode to open file with
   */
  TrajFile(const std::string &filename, int cmode) : filename(filename) {

    util::ncCheck(nc_open(filename.c_str(), cmode, &ncid), __FILE__, __LINE__);
    isOpen = true;
    bool writeable;

    // check_metadata();
    util::ncCheck(nc_inq_dimid(ncid, FRAME_NAME.c_str(), &frame_dimid),
                  __FILE__, __LINE__);
    util::ncCheck(nc_inq_dimid(ncid, NPOLYGONS_NAME.c_str(), &npolygons_dimid),
                  __FILE__, __LINE__);
    util::ncCheck(nc_inq_dimid(ncid, NVERTICES_NAME.c_str(), &nvertices_dimid),
                  __FILE__, __LINE__);
    util::ncCheck(nc_inq_dimid(ncid, SPATIAL_DIMS_NAME.c_str(), &spatial_dimid),
                  __FILE__, __LINE__);
    util::ncCheck(
        nc_inq_dimid(ncid, POLYGON_ORDER_NAME.c_str(), &polygon_order_dimid),
        __FILE__, __LINE__);

    util::ncCheck(nc_inq_varid(ncid, TOPOLOGY_DNAME.c_str(), &topology_varid),
                  __FILE__, __LINE__);
    util::ncCheck(nc_inq_varid(ncid, TIME_DNAME.c_str(), &time_varid), __FILE__,
                  __LINE__);
    util::ncCheck(nc_inq_varid(ncid, COORD_DNAME.c_str(), &coordinates_varid),
                  __FILE__, __LINE__);
    util::ncCheck(
        nc_inq_varid(ncid, MEAN_CURVATURE_DNAME.c_str(), &mean_curvature_varid),
        __FILE__, __LINE__);
    util::ncCheck(nc_inq_varid(ncid, GAUSSIAN_CURVATURE_DNAME.c_str(),
                               &gaussian_curvature_varid),
                  __FILE__, __LINE__);
  }

  /**
   * @brief Private constructor for creating a new file
   *
   * This function applies the convention. It also sets
   *
   * @param filename Path to file of interest
   * @param mesh     Mesh to store
   * @param fMode    Mode to create file with (NC_NOCLOBBER, NC_CLOBBER)
   */
  TrajFile(const std::string &filename, gcs::SurfaceMesh &mesh, int cmode)
      : filename(filename) {
    util::ncCheck(nc_create(filename.c_str(), cmode | NC_NETCDF4, &ncid),
                  __FILE__, __LINE__);

    isOpen = true;

    // Specify convention metadata
    util::ncCheck(nc_put_att_text(ncid, NC_GLOBAL, CONVENTIONS_NAME.c_str(),
                                  CONVENTIONS_VALUE.size(),
                                  CONVENTIONS_VALUE.c_str()),
                  __FILE__, __LINE__);

    util::ncCheck(nc_put_att_text(ncid, NC_GLOBAL,
                                  CONVENTIONS_VERSION_NAME.c_str(),
                                  CONVENTIONS_VERSION_VALUE.size(),
                                  CONVENTIONS_VERSION_VALUE.c_str()),
                  __FILE__, __LINE__);

    // Set dimensions
    util::ncCheck(
        nc_def_dim(ncid, FRAME_NAME.c_str(), NC_UNLIMITED, &frame_dimid),
        __FILE__, __LINE__);
    util::ncCheck(nc_def_dim(ncid, NPOLYGONS_NAME.c_str(), mesh.nFaces(),
                             &npolygons_dimid),
                  __FILE__, __LINE__);

    util::ncCheck(nc_def_dim(ncid, NVERTICES_NAME.c_str(), mesh.nVertices(),
                             &nvertices_dimid),
                  __FILE__, __LINE__);

    util::ncCheck(nc_def_dim(ncid, SPATIAL_DIMS_NAME.c_str(), SPATIAL_DIMS,
                             &spatial_dimid),
                  __FILE__, __LINE__);

    util::ncCheck(nc_def_dim(ncid, POLYGON_ORDER_NAME.c_str(), POLYGON_ORDER,
                             &polygon_order_dimid),
                  __FILE__, __LINE__);

    // Setup data blocks
    util::def_var(ncid, TOPOLOGY_DNAME, NC_UINT,
                  {npolygons_dimid, polygon_order_dimid}, &topology_varid,
                  __FILE__, __LINE__);

    util::def_var(ncid, TIME_DNAME, NC_DOUBLE, {frame_dimid}, &time_varid,
                  __FILE__, __LINE__);
    util::ncCheck(nc_put_att_text(ncid, time_varid, UNITS.c_str(),
                                  TIME_UNITS.size(), TIME_UNITS.c_str()),
                  __FILE__, __LINE__);

    util::def_var(ncid, COORD_DNAME, NC_DOUBLE,
                  {frame_dimid, nvertices_dimid, spatial_dimid},
                  &coordinates_varid, __FILE__, __LINE__);
    util::ncCheck(nc_put_att_text(ncid, coordinates_varid, UNITS.c_str(),
                                  LENGTH_UNITS.size(), LENGTH_UNITS.c_str()),
                  __FILE__, __LINE__);

    util::def_var(ncid, MEAN_CURVATURE_DNAME, NC_DOUBLE,
                  {frame_dimid, nvertices_dimid}, &mean_curvature_varid,
                  __FILE__, __LINE__);
    util::def_var(ncid, GAUSSIAN_CURVATURE_DNAME, NC_DOUBLE,
                  {frame_dimid, nvertices_dimid}, &gaussian_curvature_varid,
                  __FILE__, __LINE__);

    // Populate topology data
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
        faceMatrix = getFaceVertexMatrix(mesh);
    std::uint32_t *topodata = faceMatrix.data();
    util::ncCheck(nc_put_var(ncid, topology_varid, topodata), __FILE__,
                  __LINE__);
  }

  int ncid;

  /// Dimensions
  int frame_dimid;
  int npolygons_dimid;
  int nvertices_dimid;
  int spatial_dimid;
  int polygon_order_dimid;

  /// Data blocks
  int topology_varid;
  int time_varid;
  int coordinates_varid;
  int mean_curvature_varid;
  int gaussian_curvature_varid;

  /// Filepath to file
  std::string filename;
  /// File open status
  bool isOpen;
  /// Writeable status
  bool writeable;
};
} // namespace trajfile
} // namespace ddgsolver
#endif
