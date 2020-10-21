#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace  {

template <typename DimensionWrapper>
class InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest
 : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

using TestDimensions = ::testing::Types<bart::testing::TwoD>;
TYPED_TEST_SUITE(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
                 TestDimensions);

TYPED_TEST(InstrumentationConverterUtilityDealiiToOrderedComplexVectorTest,
       Dummy) {
  EXPECT_TRUE(false);
}

} // namespace 
