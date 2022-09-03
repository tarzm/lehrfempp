/**
 * @file
 * @brief Class representing global mesh functions in the DGFE setting
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef MESH_FUNCTION_GLOBAL_DGFE
#define MESH_FUNCTION_GLOBAL_DGFE

#include <typeinfo>

namespace lf::dgfe {

template<typename F>
class MeshFunctionGlobalDGFE {

    using F_return_type = decltype(std::declval<F>()(std::declval<Eigen::Vector2d>()));
    
public:
    MeshFunctionGlobalDGFE(const MeshFunctionGlobalDGFE&) = default;
    MeshFunctionGlobalDGFE(MeshFunctionGlobalDGFE&&) noexcept = default;
    MeshFunctionGlobalDGFE& operator=(const MeshFunctionGlobalDGFE&) = delete;
    MeshFunctionGlobalDGFE& operator=(MeshFunctionGlobalDGFE&&) = delete;

    explicit MeshFunctionGlobalDGFE(F f) : f_(std::move(f)) {}


    /**
     * @brief MeshFunction compliant evaluation operator
     * @param e the polygon in which the points lie
     * @param local the local points on which the function should be evaluated
     */
    std::vector<F_return_type> operator()(const mesh::Entity& e,
                                            const Eigen::MatrixXd& local) const {
        LF_ASSERT_MSG(e.RefEl() == lf::base::RefEl::kPolygon(),
                    "The entity must be a polygon, such that the points can be mapped from the reference bounding box into the bounding box of the polygon");
        std::vector<F_return_type> result;
        result.reserve(local.cols());
        lf::dgfe::BoundingBox box(e);
        auto global_points = box.map(local);
        for (int i = 0; i < local.cols(); i++) {
            result.push_back(f_(global_points.col(i)));
        }
        return result;                               
    }

    F lambda(){
        return f_;
    }

    ~MeshFunctionGlobalDGFE() = default;

private:
    F f_;

};


template<typename SCALAR, typename MESH_FUNC_1, typename MESH_FUNC_2> 
SCALAR L2ErrorSubTessellation(MESH_FUNC_1 g, MESH_FUNC_2 f, std::shared_ptr<const lf::mesh::Mesh> mesh_ptr, int max_degree){
    
    using return_type_1 = decltype(std::declval<MESH_FUNC_1>()(std::declval<lf::mesh::Entity>(), std::declval<Eigen::Vector2d>()));
    using return_type_2 = decltype(std::declval<MESH_FUNC_2>()(std::declval<lf::mesh::Entity>(), std::declval<Eigen::Vector2d>()));

    LF_VERIFY_MSG(typeid(return_type_1) == typeid(return_type_2), "Return types do not match");

    SCALAR error = 0.0;

    //quadrule setup
    const lf::quad::QuadRule qr_t = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), max_degree);
    // qr points
    const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
    //weights
    Eigen::VectorXd w_ref_t{qr_t.Weights()};

    for (auto cell : mesh_ptr->Entities(0)){
        //get sub-tessellation
        auto sub_tessellation = subTessellation(cell);
        //local - global mapping
        lf::dgfe::BoundingBox box(*cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){

            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

            auto dgfe_sol = g(*cell, zeta_box_t);
            auto true_sol = f(*cell, zeta_box_t);

            //sum over qr points
            for (int i = 0; i < gram_dets_t.size(); i++){
                error += (dgfe_sol[i] - true_sol[i]) * (dgfe_sol[i] - true_sol[i]) * w_ref_t[i] * gram_dets_t[i];
            }
        }
    }

    return std::sqrt(error);
}

template<typename SCALAR, typename MESH_FUNC_1, typename MESH_FUNC_2> 
SCALAR L2ErrorGradSubTessellation(MESH_FUNC_1 g, MESH_FUNC_2 f, std::shared_ptr<const lf::mesh::Mesh> mesh_ptr, int max_degree){
    
    using return_type_1 = decltype(std::declval<MESH_FUNC_1>()(std::declval<lf::mesh::Entity>(), std::declval<Eigen::Vector2d>()));
    using return_type_2 = decltype(std::declval<MESH_FUNC_2>()(std::declval<lf::mesh::Entity>(), std::declval<Eigen::Vector2d>()));

    LF_VERIFY_MSG(typeid(return_type_1) == typeid(return_type_2), "Return types do not match");
    LF_VERIFY_MSG(typeid(return_type_1) == typeid(std::vector<Eigen::Matrix<SCALAR, 2, 1>>), "Return type is not a vector");

    SCALAR error = 0.0;

    //quadrule setup
    const lf::quad::QuadRule qr_t = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), max_degree);
    // qr points
    const Eigen::MatrixXd zeta_ref_t{qr_t.Points()};
    //weights
    Eigen::VectorXd w_ref_t{qr_t.Weights()};

    for (auto cell : mesh_ptr->Entities(0)){
        //get sub-tessellation
        auto sub_tessellation = subTessellation(cell);
        //local - global mapping
        lf::dgfe::BoundingBox box(*cell);

        //loop over triangles in the sub-tessellation
        for(auto& tria_geo_ptr : sub_tessellation){

            // qr points mapped to triangle
            Eigen::MatrixXd zeta_global_t{tria_geo_ptr->Global(zeta_ref_t)};
            // qr points mapped back into reference bounding box to retrieve values
            Eigen::MatrixXd zeta_box_t{box.inverseMap(zeta_global_t)};
            //gramian determinants
            Eigen::VectorXd gram_dets_t{tria_geo_ptr->IntegrationElement(zeta_ref_t)};

            auto dgfe_sol = g(*cell, zeta_box_t);
            auto true_sol = f(*cell, zeta_box_t);

            //sum over qr points
            for (int i = 0; i < gram_dets_t.size(); i++){
                error += (dgfe_sol[i] - true_sol[i]).squaredNorm() * w_ref_t[i] * gram_dets_t[i];
            }
        }
    }

    return std::sqrt(error);
}


} //namespace lf::dgfe

#endif // MESH_FUNCTION_GLOBAL_DGFE