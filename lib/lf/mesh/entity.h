#ifndef __37e385afbd3b4b1dba8611fb71787822
#define __37e385afbd3b4b1dba8611fb71787822
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>

namespace lf::mesh {

class Entity {
public:
  /**
   * @brief The codimension of this entity w.r.t. the Mesh.dimMesh()
   * of the owning mesh manager.
   */
  virtual char Codim() const = 0;

  /**
   * @brief Return all sub entities of this entity that have the given 
   *        codimension (w.r.t. this entity!)
   * @param rel_codim The _relative co-dimension_ w.r.t. this entity
   * @return a range containing all subentities of the specified 
             _relative co-dimension_ 

   * @note For a mesh covering a manifold of dimension 2, we have the following cases
     - For a cell (co-dimension 0 entity), the cell itself is a subentity 
       of relative co-dimension 0, the edges have relative co-dimension 1, and the 
       vertices relative co-dimension 2: in this case the usual co-dimension agrees 
       with the relative co-dimension.
     - For an edge (co-dimension 1 entity), the edge itself is the only sub-entity 
       with relative co-dimension 0, and the endpoints are the sub-entitities of
       relative co-dimension 1.  
   */
  virtual base::RandomAccessRange<const Entity> SubEntities(char rel_codim) const = 0;

  /**
   * @brief Describes the geometry of this entity.
   * @return A pointer to a Geometry object that will remain valid for as long
   *         as the Mesh remains valid.
   */
  virtual geometry::Geometry* Geometry() const = 0;

  /**
   * @brief Describes the reference element type of this entity.
   * @return An object of type base::RefEl.
   */
  virtual base::RefEl RefEl() const = 0;

  /**
   * @brief Check if two entities are the same
   * @param rhs Check if this entity is the same as the rhs entity.
   * 
   * @note The behavior of this method is undefined if the rhs entity belongs 
           to a different Mesh.
   */
  virtual bool operator==(const Entity& rhs) const = 0;

  /**
   * @brief Check if two entities are different.
   */
  bool operator!=(const Entity& rhs) const { return !operator==(rhs); }

  /**
   * @brief Virtual Destructor.
   */
  virtual ~Entity() = default;
};


}

#endif // __37e385afbd3b4b1dba8611fb71787822