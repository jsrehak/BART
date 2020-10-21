#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_
#define BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_

#include <deal.II/lac/vector.h>

#include "instrumentation/converter/converter_i.h"

namespace bart::instrumentation::converter::utility {

class DealiiToOrderedComplexVector
    : public ConverterI<dealii::Vector<double>,
                        std::vector<std::complex<double>>> {
 public:
  using ComplexVector = std::vector<std::complex<double>>;
  using DealiiVector = dealii::Vector<double>;
  [[nodiscard]] ComplexVector Convert(const DealiiVector& input) const override;
};

} // namespace bart::instrumentation::converter::utility

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_
