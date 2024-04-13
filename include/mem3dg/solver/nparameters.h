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
#include <vector>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/forces.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/type_utilities.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

namespace solver {

struct nProteinParameters {
    /// Whether or not consider protein binding
    bool isProteinVariation = false;
    /// Whether conserve protein mass;
    bool isProteinConservation = false;
    /// interior point parameter for protein density
    double proteinInteriorPenalty = 0;
    /// mobility constant
    double proteinMobility = 0;
    /// binding energy per protein
    double epsilon = 0;
    /// aggregation energy constant
    double chi = 0;
    /// entropy energy constant
    double xi = 0;
    /// Smooothing coefficients
    double eta = 0;
    /// Constant of gaussian modulus vs protein density
    double Kgc = 0;
    /// Constant of bending modulus vs protein density
    double Kbc = 0;
    /// Constant of deviatoric modulus vs protein density
    double Kdc = 0;
    /// Constant of Spontaneous curvature vs protein density
    double H0c = 0;
    /// Area of protein
    double area = 0;
    /// type of relation between H0 and protein density, "linear" or "hill"
    std::string relation = "linear";
    /// precription of protein density
    std::function<EigenVectorX1d(double, EigenVectorX1d, EigenVectorX1d)>
        prescribeProteinDensityDistribution = NULL;
    /// period of updating protein density distribution
    std::size_t updateProteinDensityDistributionPeriod =
        std::numeric_limits<std::size_t>::max();
    /// prescription of center finding
    std::function<Eigen::Matrix<bool, Eigen::Dynamic, 1>(
        EigenVectorX3sr, EigenVectorX3dr, EigenVectorX1d)>
        prescribeNotableVertex = NULL;
    /// prescribeNotableVertex
    std::size_t updateNotableVertexPeriod =
        std::numeric_limits<std::size_t>::max();
};

} // namespace solver
} // namespace mem3dg