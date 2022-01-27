

#include <gtest/gtest.h>
#include <lf/base/base.h>

using namespace lf::base;

namespace lf::base::test {
TEST(RefEl, dimensionCorrect) {
  EXPECT_EQ(RefEl::kPoint().Dimension(), 0);
  EXPECT_EQ(RefEl::kSegment().Dimension(), 1);
  EXPECT_EQ(RefEl::kTria().Dimension(), 2);
  EXPECT_EQ(RefEl::kQuad().Dimension(), 2);
}

// Added ---------------------
TEST(RefEl, numberOfNodesCorrect) {
  EXPECT_EQ(RefEl::kPoint().NumNodes(), 1);
  EXPECT_EQ(RefEl::kSegment().NumNodes(), 2);
  EXPECT_EQ(RefEl::kTria().NumNodes(), 3);
  EXPECT_EQ(RefEl::kQuad().NumNodes(), 4);
}

TEST(RefEl, numSubEntitiesCorrect) {
  EXPECT_EQ(RefEl::kSegment().NumSubEntities(1), 2);  // Points for segment
  EXPECT_EQ(RefEl::kTria().NumSubEntities(1), 3);     // Segments for triangle
  EXPECT_EQ(RefEl::kTria().NumSubEntities(2), 3);     // Points for triangle
  EXPECT_EQ(RefEl::kQuad().NumSubEntities(1), 4);  // Segments for quadrilateral
  EXPECT_EQ(RefEl::kQuad().NumSubEntities(2), 4);  // Points for quadrilateral
}

TEST(RefEl, disabledPolygonCorrect) {
  EXPECT_EQ(RefEl::kPolygon().Dimension(), 2);
  EXPECT_DEATH(auto a = RefEl::kPolygon().NumNodes(), "");
  EXPECT_DEATH(auto a = RefEl::kPolygon().NodeCoords(), "");
  EXPECT_DEATH(auto a = RefEl::kPolygon().NodeCoords(), "");
  EXPECT_DEATH(auto a = RefEl::kPolygon().NumSubEntities(1), "");
  EXPECT_EQ(RefEl::kPolygon().SubType(1,1), RefEl::kSegment());
  EXPECT_EQ(RefEl::kPolygon().SubType(2,1), RefEl::kPoint());
  EXPECT_DEATH(auto a = RefEl::kPolygon().SubSubEntity2SubEntity(1, 1, 1, 1), "");
  EXPECT_EQ(RefEl::kPolygon().ToString(), "POLYGON");
}

}  // namespace lf::base::test
