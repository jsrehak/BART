#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

namespace bart::instrumentation::converter::utility {

auto DealiiToOrderedComplexVector::Convert(const DealiiVector& /*input*/) const
-> ComplexVector {
  return std::vector<std::complex<double>>();
}
} // namespace bart::instrumentation::converter::utility
