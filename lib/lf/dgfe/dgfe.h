#ifndef LF_DGFE_H
#define LF_DGFE_H

/**
 * @file
 * @brief Functionalities of Discontinuous Galerkin Finite Element Methods
 * @author Tarzis Maurer
 * @date 2022-19-04
 * @copyright ETH Zurich
 */

//#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
//#include <lf/mesh/hybrid2d/hybrid2d.h>
//#include <lf/mesh/polytopic2d/polytopic2d.h>
#include <lf/assemble/assemble.h>


#include "integration.h"
#include "bounding_box.h"
#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "dgfe_providers.h"
#include "mesh_function_dgfe.h"
#include "sub_tessellation.h"
#include "mesh_function_global.h"
#include "advection_reaction.h"
#include "discontinuity_penalization.h"
#include "advection_reaction_diffusion_rhs.h"
#include "diffusion.h"

/**
 * @brief Collects data structures and algorithms designed for discontinuous
 * galerkin finite element methods
 *
 * This namespace contains a number of classes/functions which
 * can be used to solve boundary value problems with discontinuous galerkin finite 
 * methods.
 *
 * Examples of approximation spaces that the methodsclasses in this namespace
 * can represent/handle are:
 * - Lagrangian FE of any order
 * - Hierarchical FE Spaces
 * - Approximation spaces with local p-refinement
 * - Broken spaces (e.g. for Discontinuous Galerkin Approximations)
 */
namespace lf::dgfe {


/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

// Import operators/free functions from lf::mesh::utils so we can apply them
// also to mesh functions defined in lf::fe (Argument Dependent Lookup)
using mesh::utils::operator*;
using mesh::utils::operator+;
using mesh::utils::operator-;
using mesh::utils::adjoint;
using mesh::utils::conjugate;
using mesh::utils::squaredNorm;
using mesh::utils::transpose;

}  // namespace lf::dgfe

#endif //LF_DGFE_H