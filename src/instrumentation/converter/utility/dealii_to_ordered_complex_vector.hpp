#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_
#define BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_

#include <deal.II/lac/vector.h>

#include <map>

#include "domain/definition_i.h"
#include "instrumentation/converter/converter_i.h"

namespace bart::instrumentation::converter::utility {

class DealiiToOrderedComplexVector
    : public ConverterI<dealii::Vector<double>, dealii::Vector<double>> {
 public:
  using DealiiVector = dealii::Vector<double>;
  using GlobalIndex = int;
  using OrderedIndex = int;
  using OrderingMap = std::map<GlobalIndex, OrderedIndex>;

  [[nodiscard]] DealiiVector Convert(const DealiiVector& input) const override;

  template <int dim>
  auto CalculateOrderingMap(domain::DefinitionI<dim>* domain_ptr) -> OrderingMap;

  [[nodiscard]] OrderingMap ordering_map() const noexcept { return ordering_map_; };
 private:
  OrderingMap ordering_map_{};
};

} // namespace bart::instrumentation::converter::utility

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_UTILITY_DEALII_TO_ORDERED_COMPLEX_VECTOR_HPP_
