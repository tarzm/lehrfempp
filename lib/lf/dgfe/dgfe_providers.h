/**
 * @file
 * @brief Classes supporting assembly of element matrices/vectors in the dgfe setting
 * @author Tarzis Maurer
 * @date May 22
 * @copyright ETH Zurich
*/

#ifndef DGFE_PROVIDERS_H
#define DGFE_PROVIDERS_H

#include <lf/quad/quad.h>

#include "legendre_dgfe.h"
#include "dgfe_space.h"
#include "bounding_box.h"
#include "integration.h"
#include "mesh_function_dgfe.h"
#include "mesh_function_global.h"

namespace lf::dgfe {

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/dgfe/dgfe_providers.h
 * @brief Computing the element matrix for the mass matrix using the SubTessellationIntegrator class
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref AssembleMatrixLocally().
 */
template <typename SCALAR>
class DGFEMassElementMatrixST {
    public:
        /**
         * @brief Construct a new DGFEMassElementMatrixST object
         * 
         * @param integration_max_degree the maximum degree used for the quadrature rule
         * @param legendre_max_degree maximum degree of 1D legendre polynomials: Either 1 or 2
         */
        DGFEMassElementMatrixST(unsigned max_integration_degree, unsigned max_legendre_degree) : max_integration_degree_(max_integration_degree), max_legendre_degree_(max_legendre_degree) {
            LF_VERIFY_MSG(max_legendre_degree_ == 1 || max_legendre_degree_ == 2, "Only implemented for maximum 1D legendre polynomials of degree 1 and 2");
        }

        /**
         * @brief All cells are considered active in the default implementation
         */
        [[nodiscard]] bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        /**
         * @brief main routine for the computation of element matrices
         *
         * @param cell reference to the polytopic cell for
         *        which the element matrix should be computed.
         * @return a 4x4 or a 9x9 matrix
         */
        [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> Eval(const lf::mesh::Entity &cell) const{
            unsigned matrix_size = (max_legendre_degree_ == 1) ? 4 : 9;
            Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> elem_mat(matrix_size, matrix_size);

            int i1;
            int i2;
            int j1;
            int j2;
            auto eval_lambda = [&i1, &i2, &j1, &j2](const lf::mesh::Entity &entity, Eigen::MatrixXd &local) -> std::vector<SCALAR> {
                std::vector<SCALAR> result(local.cols());
                for (int i = 0; i < local.cols(); i++){
                    result.at(i) = lf::dgfe::legendre_polynomial_2D(i1, i2, local.col(i)) * lf::dgfe::legendre_polynomial_2D(j1, j2, local.col(i));
                }
                return result;
            };

            lf::dgfe::SubTessellationIntegrator<SCALAR, decltype(eval_lambda)> integrator;
            
            //loop over trial basis funtions on cell
            for (int i = 0; i < matrix_size; i++){
                //loop over test basis functions on cell
                for (int j = 0; j < matrix_size; j++){
                    scalar_t sum = 0;
                    //definition of i1, i2, j1, j2
                    auto degrees_i = multiIndexToDegree(i, max_legendre_degree_);
                    auto degrees_j = multiIndexToDegree(j, max_legendre_degree_);
                    i1 = degrees_i.first;   //degree of x in trial basis
                    i2 = degrees_i.second;  //degree of y in trial basis
                    j1 = degrees_j.first;   //degree of x in test basis
                    j2 = degrees_j.second;  //degree of y in test basis

                    elem_mat(i,j) = integrator.integrate(cell, eval_lambda, max_integration_degree_);
                }
            }

            return elem_mat;
        }
    
    private:

        unsigned max_integration_degree_;
        unsigned max_legendre_degree_;
};


/**
 * @ingroup entity_vector_provider
 * @headerfile lf/dgfe/dgfe.h
 * @brief Local computation of general element (load) vector for scalar
 finite
 * elements; volume contributions only
 *
 * @tparam SCALAR underlying scalar type of the DGFESpace, usually double
 *
 * The underlying local linear form is
 * @f[
      v \mapsto \int_K
 f(\mathbf{x})\,\overline{v(\mathbf{x})}\,\mathrm{d}\mathbf{x}\;,
 * @f]
 * where \f$f\f$ is supposed to be a locally continuous source function.
 *
 * Computation is based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function assemble::AssembleVectorLocally().
 */
template <typename SCALAR, typename MESH_FUNCTION>
class DGFELoadElementVectorProvider {

 public:
  using ElemVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

  /** @name standard constructors
   *@{*/
  DGFELoadElementVectorProvider(const DGFELoadElementVectorProvider &) =
      delete;
  DGFELoadElementVectorProvider(DGFELoadElementVectorProvider &&) noexcept =
      default;
  DGFELoadElementVectorProvider &operator=(
      const DGFELoadElementVectorProvider &) = delete;
  DGFELoadElementVectorProvider &operator=(
      DGFELoadElementVectorProvider &&) = delete;
  /**@}*/

  /** @brief Constructor, performs precomputations
   *
   * @param fe_space specification of local shape functions
   * @param f MeshFunctionGlobalDGFE object for source function
   *
   * Uses quadrature rule of double the degree of exactness compared to the
   * degree of the finite element space.
   */
  DGFELoadElementVectorProvider(
      std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space, MESH_FUNCTION f) : f_(std::move(f)), dgfe_space_(std::move(dgfe_space)), max_legendre_degree_(dgfe_space_->MaxLegendreDegree()) {}
  /** @brief all cells are active */
  bool isActive(const lf::mesh::Entity & /*cell*/) const { return true; }

    /**
     * @brief Main method for computing the element vector
     *
     * @param cell current cell for which the element vector is desired
     * @return local load vector as column vector
     *
     */
    ElemVec Eval(const lf::mesh::Entity &cell) const {
        LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kPolygon(), "Only implemented for Polygons");

        //degree of quadrule is fixed as the maximum degree of the basis functions plus 13
        const lf::quad::QuadRule qr = qr_cache_.Get(lf::base::RefEl::kTria(), 13 + max_legendre_degree_* max_legendre_degree_);

        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);
        // qr points
        const Eigen::MatrixXd zeta_ref{qr.Points()};
        //weights
        Eigen::VectorXd w_ref{qr.Weights()};

        //initialize element vector
        unsigned vector_size = (max_legendre_degree_ == 1) ? 4 : 9;
        Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> elem_vec(vector_size, 1);
        elem_vec.setZero();

        lf::dgfe::BoundingBox box(cell);

        int i1;
        int i2;
        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global{tria_geo_ptr->Global(zeta_ref)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box{box.inverseMap(zeta_global)};
            //gramian determinants
            Eigen::VectorXd gram_dets{tria_geo_ptr->IntegrationElement(zeta_ref)};
            //loop over basis functions
            for (int basis = 0; basis < vector_size; basis++){

                auto basis_degrees = multiIndexToDegree(basis, max_legendre_degree_);
                i1 = basis_degrees.first;
                i2 = basis_degrees.second;
                //sum over qr points
                for (int i = 0; i < qr.Points().cols(); i++){

                    //Note: functor calls entity as argument
                    elem_vec[basis] += w_ref[i] * f_(cell, zeta_box.col(i))[0] * legendre_polynomial_2D(i1, i2, zeta_box.col(i)) * gram_dets[i];
                }
            }
        }
        return elem_vec;
    }

  ~DGFELoadElementVectorProvider() = default;

 private:
  /** @brief An object providing the source function */
  MESH_FUNCTION f_;

  std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_;

  lf::quad::QuadRuleCache qr_cache_;

  unsigned max_legendre_degree_;
};

template <typename SCALAR, typename MESH_FUNCTION>
class L2ProjectionSqrtANablaBasisLoadVector {

 public:

  /** @name standard constructors
   *@{*/
  L2ProjectionSqrtANablaBasisLoadVector(const L2ProjectionSqrtANablaBasisLoadVector &) =
      delete;
  L2ProjectionSqrtANablaBasisLoadVector(L2ProjectionSqrtANablaBasisLoadVector &&) noexcept =
      default;
  L2ProjectionSqrtANablaBasisLoadVector &operator=(
      const L2ProjectionSqrtANablaBasisLoadVector &) = delete;
  L2ProjectionSqrtANablaBasisLoadVector &operator=(
      L2ProjectionSqrtANablaBasisLoadVector &&) = delete;
  /**@}*/

  /** @brief Constructor, performs precomputations
   *
   * @param fe_space specification of local shape functions
   * @param f MeshFunctionGlobalDGFE object for source function
   *
   * Uses quadrature rule of double the degree of exactness compared to the
   * degree of the finite element space.
   */
  L2ProjectionSqrtANablaBasisLoadVector(
      std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space, MESH_FUNCTION a_coeff, unsigned dim, unsigned basis, unsigned max_integration_degree) : a_coeff_(a_coeff), dgfe_space_(std::move(dgfe_space)),
                             max_legendre_degree_(dgfe_space_->MaxLegendreDegree()), dim_(dim), basis_(basis), max_integration_degree_(max_integration_degree) {
                                LF_VERIFY_MSG(dim == 1 || dim == 0, "Desired dimension of gradient has to be either 0 or 1");
                                LF_VERIFY_MSG(0 <= basis_ && basis <= (max_legendre_degree_ + 1) * (max_legendre_degree_ + 1),
                                                 "basis function index needs to be within the range of number of local dofs");
                             }

  /** @brief all cells are active */
  bool isActive(const lf::mesh::Entity & /*cell*/) const { return true; }

    /**
     * @brief Main method for computing the element vector
     *
     * @param cell current cell for which the element vector is desired
     * @return local load vector as column vector 
     *
     */
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> Eval(const lf::mesh::Entity &cell) const {
        LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kPolygon(), "Only implemented for Polygons");

        //degree of quadrule is fixed as the maximum degree of the basis functions plus 13
        const lf::quad::QuadRule qr = qr_cache_.Get(lf::base::RefEl::kTria(), max_integration_degree_);

        //get sub-tessellation
        auto sub_tessellation = subTessellation(&cell);
        // qr points
        const Eigen::MatrixXd zeta_ref{qr.Points()};
        //weights
        Eigen::VectorXd w_ref{qr.Weights()};

        //initialize element vector
        unsigned vector_size = (max_legendre_degree_ == 1) ? 4 : 9;
        Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> elem_vec(vector_size, 1);
        elem_vec.setZero();

        lf::dgfe::BoundingBox box(cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){
            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global{tria_geo_ptr->Global(zeta_ref)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box{box.inverseMap(zeta_global)};
            //gramian determinants
            Eigen::VectorXd gram_dets{tria_geo_ptr->IntegrationElement(zeta_ref)};
            
            //diffusion tensor evaluated at qr points
            auto a_eval = a_coeff_(cell, zeta_box);

        
            //loop over basis functions
            for (int basis_test = 0; basis_test < vector_size; basis_test++){

                //sum over qr points
                for (int i = 0; i < qr.Points().cols(); i++){

                    Eigen::Vector2d nabla_basis{legendre_basis_dx(basis_, max_legendre_degree_, zeta_box.col(i)) * box.inverseJacobi(0),
                                                legendre_basis_dy(basis_, max_legendre_degree_, zeta_box.col(i)) * box.inverseJacobi(1)};
                
                    elem_vec(basis_test) += ((a_eval[i].sqrt()).row(dim_)).dot(nabla_basis) * legendre_basis(basis_test, max_legendre_degree_, zeta_box.col(i)) * w_ref[i] * gram_dets[i];
                }
            }
        }
        return elem_vec;
    }

 private:
    /** @brief An object providing the diffusion tensor */
    MESH_FUNCTION a_coeff_;

    std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_;

    lf::quad::QuadRuleCache qr_cache_;

    unsigned max_integration_degree_;

    unsigned max_legendre_degree_;
    //spefifying whether 0th or 1st dimension of the gradient vector is calculated
    unsigned dim_;
    //index of the basis function whose gradient should be projected
    unsigned basis_;
};



} //namespace lf::dgfe


#endif //DGFE_PROVIDERS_H