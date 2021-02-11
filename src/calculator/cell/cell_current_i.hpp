#ifndef BART_SRC_CALCULATOR_CELL_CELL_CURRENT_I_HPP_
#define BART_SRC_CALCULATOR_CELL_CELL_CURRENT_I_HPP_

#include <memory>

#include <deal.II/lac/vector.h>

#include "domain/domain_types.h"
#include "quadrature/quadrature_types.h"
#include "utility/has_description.h"
#include "utility/has_dependencies.h"

namespace bart::calculator::cell {

template <int dim>
class CellCurrentI : public utility::HasDescription, public utility::HasDependencies{
 public:
  using QuadraturePointIndex = quadrature::QuadraturePointIndex;
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<Vector>;
  using AngularFluxMap = std::map<QuadraturePointIndex, VectorPtr>;

  using CellPtr = domain::CellPtr<dim>;
  virtual ~CellCurrentI() = default;
  virtual auto CurrentAtQuadrature(const CellPtr&, const AngularFluxMap&) const -> std::vector<Vector> = 0;
  virtual auto EddingtonCurrentAtQuadrature(const CellPtr&, const AngularFluxMap&, const double sigma_t) const
  -> std::vector<Vector> = 0;
};

} // namespace bart::calculator::cell

#endif //BART_SRC_CALCULATOR_CELL_CELL_CURRENT_I_HPP_
