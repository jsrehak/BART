#include "calculator/cell/cell_current.hpp"

namespace bart::calculator::cell {

template<int dim>
CellCurrent<dim>::CellCurrent(std::shared_ptr<FiniteElement> finite_element_ptr,
                              std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr)
    : finite_element_ptr_(finite_element_ptr), angular_flux_integrator_ptr_(angular_flux_integrator_ptr) {
  std::string call_location{"CellCurrent constructor"};
  this->AssertPointerNotNull(finite_element_ptr_.get(), "finite element", call_location);
  this->AssertPointerNotNull(angular_flux_integrator_ptr_.get(), "angular flux integrator", call_location);
}

template<int dim>
auto CellCurrent<dim>::CurrentAtQuadrature(const CellPtr& /*cell_ptr*/,
                                           const AngularFluxMap& /*angular_flux_map*/) const -> std::vector<Vector> {
  return std::vector<Vector>();
}

template<int dim>
auto CellCurrent<dim>::EddingtonCurrentAtQuadrature(const CellPtr& cell_ptr,
                                                    const AngularFluxMap& /*angular_flux_map*/,
                                                    const double /*sigma_t*/) const -> std::vector<Vector> {
  return std::vector<Vector>();
}

template class CellCurrent<1>;
template class CellCurrent<2>;
template class CellCurrent<3>;

} // namespace bart::calculator::cell
