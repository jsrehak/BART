#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

namespace bart::instrumentation::converter::utility {

auto DealiiToOrderedComplexVector::Convert(const DealiiVector& /*input*/) const
-> ComplexVector {
  return std::vector<std::complex<double>>();
}

template<int dim>
DealiiToOrderedComplexVector::OrderingMap DealiiToOrderedComplexVector::CalculateOrderingMap(
    domain::DefinitionI<dim> *domain_ptr) {
  return bart::instrumentation::converter::utility::DealiiToOrderedComplexVector::OrderingMap();
}

} // namespace bart::instrumentation::converter::utility
