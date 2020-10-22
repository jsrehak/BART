#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

#include <deal.II/base/geometry_info.h>

#include "domain/tests/definition_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace  {



template <typename DimensionWrapper>
class InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest
 : public ::testing::Test,
   public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainType = bart::domain::DefinitionMock<dim>;
  using TestConverterType = bart::instrumentation::converter::utility::DealiiToOrderedComplexVector;

  InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest()
      : bart::testing::DealiiTestDomain<DimensionWrapper::value>(1.0, 4) {}
  std::shared_ptr<DomainType> domain_ptr_;
  void SetUp() override;
};

template<typename DimensionWrapper>
void InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest<
    DimensionWrapper>::SetUp() {
  domain_ptr_ = std::make_shared<DomainType>();
  this->SetUpDealii();
}

using TestDimensions = ::testing::Types<bart::testing::TwoD>;
TYPED_TEST_SUITE(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
                 TestDimensions);

::testing::AssertionResult IsOrdered(std::vector<dealii::Point<2>> points) {
  for (std::size_t i = 1; i < points.size(); ++i) {
    auto& last_point{ points.at(i - 1) };
    auto& this_point{ points.at(i) };
    double this_y = this_point[1], this_x = this_point[0],
    last_y = last_point[1], last_x = last_point[0];

    std::cout << "Point " << i << " y = " << this_y << " x = " << this_x << "\n";

    if (this_y < last_y) {
      // If y-value is less than previous, immppediate failure
      return ::testing::AssertionFailure()
          << "Index: " << i << " has y value " << this_y
          << " but previous point has y value " << last_y;
    } else if (this_y == last_y) {
      // If y-value is equal, check x value;
      if (this_x < last_x) {
        return ::testing::AssertionFailure()
            << "Index: " << i << " has x value " << this_x
            << " but previous point has x value " << last_x;
      }
    }
  }
  return ::testing::AssertionSuccess();
}

TYPED_TEST(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
       Dummy) {
  constexpr int dim = this->dim;
  using Point = dealii::Point<dim>;

  std::vector<Point> points(this->locally_owned_dofs_.size());

  constexpr int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
  for (const auto& cell : this->cells_) {
    for (int index = 0; index < vertices_per_cell; ++index) {
      points.at(cell->vertex_index(index)) = cell->vertex(index);
    }
  }

  EXPECT_TRUE(IsOrdered(points));
}

} // namespace 
