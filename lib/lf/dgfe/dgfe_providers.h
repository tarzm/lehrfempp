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

//polyomial expansion of a function in 2D
//format: {[coefficient, (degree x, degree y)], [coefficient, (degree x, degree y)], ... }
using Polynomial = std::vector<std::pair<scalar_t, std::pair<size_type, size_type>>>;


/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/dgfe/dgfe_providers.h
 * @brief Computing the element matrix for the mass matrix
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref AssembleMatrixLocally().
 */
class DGFEO1MassElementMatrix {
    public:
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
         * @return a 4x4 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, 4, 4> Eval(const lf::mesh::Entity &cell) const;
};

class DGFEO1LocalLoadVector {
    public:
        using ElemVec = Eigen::Matrix<scalar_t, 4, 1>;

        DGFEO1LocalLoadVector(Polynomial polynomial) : polynomial_(polynomial) {}

        bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        ElemVec Eval(const lf::mesh::Entity &cell) const;

    private:
        Polynomial polynomial_;


};

class DGFEO2MassElementMatrix {
    public:
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
         * @return a 9x9 matrix
         */
        [[nodiscard]] Eigen::Matrix<scalar_t, 9, 9> Eval(const lf::mesh::Entity &cell) const;
};

class DGFEO2LocalLoadVector {
    public:
        using ElemVec = Eigen::Matrix<scalar_t, 9, 1>;
        

        DGFEO2LocalLoadVector(Polynomial polynomial) : polynomial_(polynomial) {}

        bool isActive(const lf::mesh::Entity & /*cell*/) const {
            return true;
        }

        ElemVec Eval(const lf::mesh::Entity &cell) const;

    private:
        Polynomial polynomial_;


};

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
template <typename SCALAR, typename FUNCTOR>
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
      std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space, lf::dgfe::MeshFunctionGlobalDGFE<FUNCTOR> f) : f_(std::move(f)), dgfe_space_(std::move(dgfe_space)), max_legendre_degree_(dgfe_space_->MaxLegendreDegree()) {}
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
        const lf::quad::QuadRule qr = qr_cache_.Get(lf::base::RefEl::kTria(), 12 + max_legendre_degree_* max_legendre_degree_);

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
  lf::dgfe::MeshFunctionGlobalDGFE<FUNCTOR> f_;

  std::shared_ptr<const lf::dgfe::DGFESpace> dgfe_space_;

  lf::quad::QuadRuleCache qr_cache_;

  unsigned max_legendre_degree_;
};




} //namespace lf::dgfe


#endif //DGFE_PROVIDERS_H