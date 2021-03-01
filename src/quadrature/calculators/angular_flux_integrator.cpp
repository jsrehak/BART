#include "quadrature/calculators/angular_flux_integrator.hpp"
#include "quadrature/calculators/quadrature_calculators_factories.hpp"

namespace bart::quadrature::calculators {

template<int dim>
AngularFluxIntegrator<dim>::AngularFluxIntegrator(std::shared_ptr<QuadratureSet> quadrature_set_ptr)
    : quadrature_set_ptr_(quadrature_set_ptr) {
  this->AssertPointerNotNull(quadrature_set_ptr_.get(), "quadrature set", "DriftDiffusionIntegratedFlux constructor");
}

template<int dim>
bool AngularFluxIntegrator<dim>::is_registered_ =
    AngularFluxIntegrator<dim>::Factory::get().RegisterConstructor(
        AngularFluxIntegratorName::kDefaultImplementation,
        [](std::shared_ptr<QuadratureSet> quadrature_set_ptr) -> std::unique_ptr<AngularFluxIntegratorI> {
          return std::make_unique<AngularFluxIntegrator<dim>>(quadrature_set_ptr); });

template<int dim>
auto AngularFluxIntegrator<dim>::EddingtonCurrent(const GradientMap& gradient_map,
                                                  const double sigma_t) const -> std::vector<Vector>{
  using AngleIndex = quadrature::QuadraturePointIndex;
  const auto max_index{ gradient_map.cbegin()->second.size() };
  const GradientMap::size_type n_quadrature_points{ this->quadrature_set_ptr_->size() };
  std::vector<Vector> return_vector(max_index);

  for (GradientMap::size_type angle_index = 0; angle_index < n_quadrature_points; ++angle_index) {
    auto quadrature_point = quadrature_set_ptr_->GetQuadraturePoint(AngleIndex(angle_index));
    const double weight{ quadrature_point->weight() };
    const auto position{ quadrature_point->cartesian_position_tensor() };
    const double sum_constant{ weight * position * position /sigma_t };
    for (std::vector<Vector>::size_type index = 0; index < max_index; ++index) {
      auto& index_vector = return_vector.at(index);
      if (index_vector.size() == 0)
        index_vector.reinit(dim);
      index_vector.add(sum_constant, gradient_map.at(AngleIndex(angle_index)).at(index));
    }
  }
  return return_vector;
}

template<int dim>
auto AngularFluxIntegrator<dim>::NetCurrent(const VectorMap& angular_flux_map) const -> std::vector<Vector> {
  const auto n_dofs{ angular_flux_map.cbegin()->second->size() };
  std::vector<Vector> return_vector;
  for (unsigned int i = 0; i < n_dofs; ++i) {
    return_vector.push_back(NetCurrent(angular_flux_map, DegreeOfFreedom(i)));
  }
  return return_vector;
}

template<int dim>
auto AngularFluxIntegrator<dim>::NetCurrent(const VectorMap& angular_flux_map,
                                            const DegreeOfFreedom degree_of_freedom) const -> Vector {
  using Index = quadrature::QuadraturePointIndex;
  const VectorMap::size_type n_quadrature_points{ this->quadrature_set_ptr_->size() };

  Vector result_vector(dim);

  for (VectorMap::size_type i = 0; i < n_quadrature_points; ++i) {
    auto& quadrature_point = *quadrature_set_ptr_->GetQuadraturePoint(Index(i));
    const double weight{ quadrature_point.weight() };
    const auto position{ quadrature_point.cartesian_position_tensor() };
    Vector position_vector(position.begin_raw(), position.end_raw());
    const auto angular_flux_at_dof = (*angular_flux_map.at(Index(i)))[degree_of_freedom.get()];
    result_vector.add(weight * angular_flux_at_dof, position_vector);
  }

  return result_vector;
}

template<int dim>
auto AngularFluxIntegrator<dim>::DirectionalCurrent(const VectorMap& angular_flux_map,
                                                    const Vector normal_vector) const -> std::vector<double> {
  const auto n_dofs{ angular_flux_map.cbegin()->second->size() };
  std::vector<double> return_vector;
  for (unsigned int i = 0; i < n_dofs; ++i) {
    return_vector.push_back(DirectionalCurrent(angular_flux_map, normal_vector, DegreeOfFreedom(i)));
  }
  return return_vector;
}

template<int dim>
auto AngularFluxIntegrator<dim>::DirectionalCurrent(const VectorMap& angular_flux_map,
                                                    const Vector normal,
                                                    const DegreeOfFreedom degree_of_freedom) const -> double {
  using Index = quadrature::QuadraturePointIndex;
  const VectorMap::size_type n_quadrature_points{ this->quadrature_set_ptr_->size() };

  double result{ 0 };

  for (VectorMap::size_type i = 0; i < n_quadrature_points; ++i) {
    auto& quadrature_point = *quadrature_set_ptr_->GetQuadraturePoint(Index(i));
    const double weight{ quadrature_point.weight() };
    const auto position{ quadrature_point.cartesian_position_tensor() };
    Vector position_vector(position.begin_raw(), position.end_raw());
    AssertThrow(position_vector.size() == normal.size(), dealii::ExcMessage("Error in DirectionalCurrentIntegration, "
                                                                            "normal vector is incorrect size"))
    const double omega_dot_normal = position_vector * normal;
    if (omega_dot_normal >= 0) {
      const auto angular_flux_at_dof = (*angular_flux_map.at(Index(i)))[degree_of_freedom.get()];
      result += weight * omega_dot_normal * angular_flux_at_dof;
    }
  }

  return result;
}

template<int dim>
auto AngularFluxIntegrator<dim>::DirectionalFlux(const VectorMap& angular_flux_map,
                                                 const Vector normal_vector) const -> std::vector<double> {
  const auto n_dofs{ angular_flux_map.cbegin()->second->size() };
  std::vector<double> return_vector;
  for (unsigned int i = 0; i < n_dofs; ++i) {
    return_vector.push_back(DirectionalFlux(angular_flux_map, normal_vector, DegreeOfFreedom(i)));
  }
  return return_vector;
}

template<int dim>
auto AngularFluxIntegrator<dim>::DirectionalFlux(const VectorMap& angular_flux_map,
                                                 const Vector normal,
                                                 DegreeOfFreedom degree_of_freedom) const -> double {
  using Index = quadrature::QuadraturePointIndex;
  const VectorMap::size_type n_quadrature_points{ this->quadrature_set_ptr_->size() };

  double result{ 0 };

  for (VectorMap::size_type i = 0; i < n_quadrature_points; ++i) {
    auto& quadrature_point = *quadrature_set_ptr_->GetQuadraturePoint(Index(i));
    const double weight{ quadrature_point.weight() };
    const auto position{ quadrature_point.cartesian_position_tensor() };
    Vector position_vector(position.begin_raw(), position.end_raw());
    AssertThrow(position_vector.size() == normal.size(), dealii::ExcMessage("Error in DirectionalFluxIntegration, "
                                                                            "normal vector is incorrect size"))
    const double omega_dot_normal = position_vector * normal;
    if (omega_dot_normal >= 0) {
      const auto angular_flux_at_dof = (*angular_flux_map.at(Index(i)))[degree_of_freedom.get()];
      result += weight * angular_flux_at_dof;
    }
  }

  return result;
}

template class AngularFluxIntegrator<1>;
template class AngularFluxIntegrator<2>;
template class AngularFluxIntegrator<3>;

} // namespace bart::quadrature::calculators
