#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

#include <deal.II/base/geometry_info.h>

#include "domain/tests/definition_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using ::testing::Return, ::testing::ContainerEq;

namespace test_helpers = bart::test_helpers;

template <typename DimensionWrapper>
class InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest
 : public ::testing::Test,
   public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainType = bart::domain::DefinitionMock<dim>;
  using TestConverterType = bart::instrumentation::converter::utility::DealiiToOrderedComplexVector;

  InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest()
      : bart::testing::DealiiTestDomain<DimensionWrapper::value>(1.0, 3) {}
  std::shared_ptr<DomainType> domain_ptr_;
  void SetUp() override;
};

template<typename DimensionWrapper>
void InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest<
    DimensionWrapper>::SetUp() {
  domain_ptr_ = std::make_shared<DomainType>();
  this->SetUpDealii();
}

using TestDimensions = ::testing::Types<bart::testing::TwoD, bart::testing::ThreeD>;
TYPED_TEST_SUITE(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
                 TestDimensions);

template <int dim>
::testing::AssertionResult IsOrdered(std::vector<dealii::Point<dim>> points) {
  for (std::size_t i = 1; i < points.size(); ++i) {
    auto &last_point{points.at(i - 1)};
    auto &this_point{points.at(i)};

    if (dim == 3) {
      double this_z = this_point[2], last_z = last_point[2];
      if (this_z < last_z) {
        return ::testing::AssertionFailure()
            << "Index: " << i << " has z value " << this_z
            << " but previous point has z value " << last_z;
      } else if (this_z > last_z) {
        break;
      }
    }
    if (dim >= 2) {
      double this_y = this_point[1], last_y = last_point[1];
      if (this_y < last_y) {
        // If y-value is less than previous, fail
        return ::testing::AssertionFailure()
            << "Index: " << i << " has y value " << this_y
            << " but previous point has y value " << last_y;
      } else if (this_y > last_y) {
        break;
      }
    }
    // If x-value is less than previous, fail
    double this_x = this_point[0], last_x = last_point[0];
    if (this_x < last_x) {
      return ::testing::AssertionFailure()
          << "Index: " << i << " has x value " << this_x
          << " but previous point has x value " << last_x;
    }
  }
  return ::testing::AssertionSuccess();
}

TYPED_TEST(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
           OrderingMap) {
  constexpr int dim = this->dim;
  using Point = dealii::Point<dim>;
  using TestConverterType = typename InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest<TypeParam>::TestConverterType;

  std::vector<Point> points(this->locally_owned_dofs_.size());
  constexpr int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
  for (const auto& cell : this->cells_) {
    for (int index = 0; index < vertices_per_cell; ++index) {
      points.at(cell->vertex_dof_index(index, 0)) = cell->vertex(index);
    }
  }
  TestConverterType test_converter;

  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(Return(this->cells_));
  auto ordering_map = test_converter.CalculateOrderingMap(this->domain_ptr_.get());

  ASSERT_EQ(ordering_map.size(), points.size());
  std::vector<Point> ordered_points(this->locally_owned_dofs_.size());
  for (const auto& [global_index, ordered_index] : ordering_map) {
    ordered_points.at(ordered_index) = points.at(global_index);
  }
  EXPECT_TRUE(IsOrdered(ordered_points));
  EXPECT_THAT(ordering_map, ContainerEq(test_converter.ordering_map()));
}

TYPED_TEST(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
           ConvertNoOrderingMap) {
  using DealiiVector = dealii::Vector<double>;
  using TestConverterType = typename InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest<TypeParam>::TestConverterType;
  TestConverterType test_converter;
  EXPECT_ANY_THROW({
    test_converter.Convert(DealiiVector{});
  });
}

TYPED_TEST(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
           Convert) {
  using DealiiVector = dealii::Vector<double>;
  using TestConverterType = typename InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest<TypeParam>::TestConverterType;

  DealiiVector to_convert(this->locally_owned_dofs_.size());
  for (auto& value : to_convert)
    value = test_helpers::RandomDouble(-100, 100);

  TestConverterType test_converter;
  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(Return(this->cells_));
  auto ordering_map = test_converter.CalculateOrderingMap(this->domain_ptr_.get());

  DealiiVector ordered_vector(this->locally_owned_dofs_.size());
  for (auto& [global_index, ordered_index] : ordering_map) {
    ordered_vector[ordered_index] = to_convert[global_index];
  }

  auto converted_vector = test_converter.Convert(to_convert);
  ASSERT_NE(converted_vector.size(), 0);
  EXPECT_EQ(converted_vector, ordered_vector);
}

} // namespace 
