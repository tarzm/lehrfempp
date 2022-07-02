/**
 * @file
 * @brief Class representing global mesh functions in the DGFE setting
 * @author Tarzis Maurer
 * @date June 22
 * @copyright ETH Zurich
*/

#ifndef MESH_FUNCTION_GLOBAL_DGFE
#define MESH_FUNCTION_GLOBAL_DGFE

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

    ~MeshFunctionGlobalDGFE() = default;

private:
    F f_;

};

} //namespace lf::dgfe

#endif // MESH_FUNCTION_GLOBAL_DGFE