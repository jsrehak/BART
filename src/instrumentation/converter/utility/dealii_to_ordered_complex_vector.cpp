#include "instrumentation/converter/utility/dealii_to_ordered_complex_vector.hpp"

#include <algorithm>
#include <utility>

#include <deal.II/base/geometry_info.h>

namespace bart::instrumentation::converter::utility {

auto DealiiToOrderedComplexVector::Convert(const DealiiVector& input) const
-> DealiiVector {
  AssertThrow(ordering_map_.size() > 0,
      dealii::ExcMessage("error in DealiiToOrderedComplexVector::Convert, ordering_map is empty, have you called "
                         "CalculateOrderingMap?"))
  AssertThrow(ordering_map_.size() == input.size(),
      dealii::ExcMessage("Error in DealiiToOrderedComplexVector::Convert, ordering_map is not the same size as "
                         " input."))
  DealiiVector ordered_vector(input.size());
  for (auto& [global_index, ordered_index] : ordering_map_) {
    ordered_vector[ordered_index] = input[global_index];
  }

  return ordered_vector;
}

template<int dim>
auto DealiiToOrderedComplexVector::CalculateOrderingMap(domain::DefinitionI<dim> *domain_ptr) -> OrderingMap {
  using Point = dealii::Point<dim>;
  constexpr int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
  auto cell_global_index_and_point = [](domain::CellPtr<dim> cell, const int index) -> std::pair<Point, int> {
    return std::pair{cell->vertex(index), cell->vertex_dof_index(index, 0)};
  };
  struct PointCompare {
    bool operator()(const Point &lhs, const Point &rhs) const {
      // Returns true if lhs < rhs
      if (dim == 3) {
        double lhs_z = lhs[2], rhs_z = rhs[2];
        if (lhs_z != rhs_z)
          return lhs_z < rhs_z;
      }
      if (dim >= 2) {
        double lhs_y = lhs[1], rhs_y = rhs[1];
        if (lhs_y != rhs_y)
          return lhs_y < rhs_y;
      }
      double lhs_x = lhs[0], rhs_x = rhs[0];
      return lhs_x < rhs_x;
    };
  };

  std::map<Point, int, PointCompare> global_index_and_point_pairs;

  for (const auto& cell : domain_ptr->Cells()) {
    for (int index = 0; index < vertices_per_cell; ++index) {
      global_index_and_point_pairs.emplace(cell_global_index_and_point(cell, index));
    }
  }

  for (std::size_t ordered_index{ 0 }; auto global_index_and_point_pair : global_index_and_point_pairs) {
    ordering_map_.emplace(global_index_and_point_pair.second, ordered_index);
    ++ordered_index;
  }

  return ordering_map_;
}

template auto DealiiToOrderedComplexVector::CalculateOrderingMap(domain::DefinitionI<2> *domain_ptr) -> OrderingMap;
template auto DealiiToOrderedComplexVector::CalculateOrderingMap(domain::DefinitionI<3> *domain_ptr) -> OrderingMap;


} // namespace bart::instrumentation::converter::utility
