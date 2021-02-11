#ifndef BART_SRC_CALCULATOR_CELL_CELL_CURRENT_HPP_
#define BART_SRC_CALCULATOR_CELL_CELL_CURRENT_HPP_

#include "domain/finite_element/finite_element_i.hpp"
#include "quadrature/calculators/angular_flux_integrator_i.hpp"
#include "calculator/cell/cell_current_i.hpp"

namespace bart::calculator::cell {

template <int dim>
class CellCurrent : public CellCurrentI<dim> {
 public:
  using typename CellCurrentI<dim>::CellPtr, typename CellCurrentI<dim>::AngularFluxMap;
  using typename CellCurrentI<dim>::VectorPtr, typename CellCurrentI<dim>::QuadraturePointIndex;
  using typename CellCurrentI<dim>::Vector;

  using FiniteElement = domain::finite_element::FiniteElementI<dim>;
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorI;

  CellCurrent(std::shared_ptr<FiniteElement>, std::shared_ptr<AngularFluxIntegrator>);

  [[nodiscard]] auto CurrentAtQuadrature(const CellPtr&, const AngularFluxMap&) const -> std::vector<Vector> override;
  [[nodiscard]] auto EddingtonCurrentAtQuadrature(const CellPtr&, const AngularFluxMap&, const double sigma_t) const
  -> std::vector<Vector> override;

  auto finite_element_ptr() -> FiniteElement* { return finite_element_ptr_.get(); }
  auto angular_flux_integrator_ptr() -> AngularFluxIntegrator* { return angular_flux_integrator_ptr_.get(); }
 private:
  std::shared_ptr<FiniteElement> finite_element_ptr_{ nullptr };
  std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr_{ nullptr };
};

} // namespace bart::calculator::cell

#endif //BART_SRC_CALCULATOR_CELL_CELL_CURRENT_HPP_
