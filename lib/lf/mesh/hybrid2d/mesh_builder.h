#ifndef __e98a803fac5b430a8ff634ceb2f809a1
#define __e98a803fac5b430a8ff634ceb2f809a1

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements mesh::MeshBuilder interface and can be used to construct
 *        a hybrid mesh with `dimMesh=2`.
 *
 * A planar triangular mesh with affine triangles as cells can be completely specified
 * by giving the list of vertex coordinates and a list of triangles in the form of a 
 * 3-tuple of vertex indices. 
 */
class MeshBuilder : public mesh::MeshBuilder {
  
public:

  /**
   * @brief Construct a new builder that can be used to construct a new hybrid2d
   *        mesh.
   * @param dim_world The dimension of the euclidean space in which the
   *                  mesh is embedded.
   */
  MeshBuilder(dim_t dim_world) : dim_world_(dim_world), built_(false) {}

  /** @copydoc Mesh::DimWorld */
  dim_t DimWorld() const override {return dim_world_; }

  /** 
   * @brief 2D hybrid meshes are meant to model 2D manifolds
   *
   */
  dim_t DimMesh() const override {return 2; }

  /** 
   * @brief register the coordinates of another point
   * @param coord a dynamic Eigen vector of size DimWorld containing node coordinates
   *
   */
  size_type AddPoint(coord_t coord) override;

  /* 
   * @brief register another triangle
   * @param nodes a 3-tuple of node *indices*
   * @param geometry a description of the shape of a cell
   *
   */
  size_type AddElement(const base::ForwardRange<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) override;

  /** 
   * @brief actual construction of the mesh 
   *
   */
  std::unique_ptr<mesh::Mesh> Build() override;
private:
  dim_t dim_world_; // dimension of ambient space
  bool built_;      
  std::vector<Eigen::VectorXd> nodes_;
  std::vector<std::tuple<std::vector<size_type>,
    std::unique_ptr<geometry::Geometry>>> elements_;
};

/**
 * @brief Implements a `MeshBuilder` that generates a triangular structured mesh
 *
 * Mesh generator for a triangular tensor product mesh covering a rectangle.
 *
 * A triangular tensor product mesh of a rectangular domain is built by
 * subdividing the domain into equal squares and splitting each of them
 * along the diagonal.
 *
 * The nodes and cells of the resulting mesh are numbered lexikographically
 * from bottom left to  top right. The super-diagonal cells have even indices
 * the sub-diagonal cells odd indices.
 * 
 * The horizontal edges are numbered first, the vertical edges next; both groups
 * lexikographically.
 */
class TPTriagMeshBuilder: public mesh::MeshBuilder {
public:
  /**
   * @brief Consrtructor
   *
   */
  TPTriagMeshBuilder():no_of_x_cells_(0),no_of_y_cells_(0) {} 
   
  /** @copydoc Mesh::DimWorld */
  dim_t DimWorld() const override {return 2; }

  /** 
   * @brief 2D hybrid meshes are meant to model 2D manifolds
   *
   */
  dim_t DimMesh() const override {return 2; }

  /**
   * @brief Initialization methods
   *
   */
  template <typename VECTOR>
  TPTriagMeshBuilder &setBottomLeftCorner(const VECTOR &&blc) {
    bottom_left_corner_ << blc[0],blc[1]; return *this;
  }
  template <typename VECTOR>
  TPTriagMeshBuilder &setTopRightCorner(const VECTOR &&trc) {
    bottom_left_corner_ << trc[0],trc[1]; return *this;
  }
  TPTriagMeshBuilder &setNoXCells(size_type nxc) { 
    no_of_x_cells_ = nxc; return *this;
  }
  TPTriagMeshBuilder &setNoYCells(size_type nyc) { 
    no_of_x_cells_ = nyc; return *this;
  }
  
  /** 
   * @brief actual construction of the mesh 
   *
   */
  std::unique_ptr<mesh::Mesh> Build() override;
private:
  /**
   * @brief vertex index from grid position
   */
  inline size_type VertexIndex(size_type i,size_type j) const {
    return (i+j*(no_of_x_cells_+1));
  }
  
private:
  Eigen::Vector2d bottom_left_corner_,top_right_corner_;
  size_type no_of_x_cells_,no_of_y_cells_;
  
};

}



#endif // __e98a803fac5b430a8ff634ceb2f809a1
