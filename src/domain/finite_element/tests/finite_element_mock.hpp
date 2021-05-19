#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_HPP_
#define BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_HPP_

#include "domain/finite_element/finite_element_i.hpp"

#include <deal.II/fe/fe_values.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart::domain::finite_element {

template <int dim>
class FiniteElementMock : public FiniteElementI<dim> {
 public:
  using typename FiniteElementI<dim>::DealiiVector, typename FiniteElementI<dim>::Tensor;
  MOCK_METHOD(int, polynomial_degree, (), (const, override));
  MOCK_METHOD(int, dofs_per_cell, (), (const, override));
  MOCK_METHOD(int, n_cell_quad_pts, (), (const, override));
  MOCK_METHOD(int, n_face_quad_pts, (), (const, override));
  MOCK_METHOD(bool, SetCell, (const domain::CellPtr<dim> &));
  MOCK_METHOD(bool, SetFace, (const domain::CellPtr<dim> &to_set, const domain::FaceIndex), (override));
  MOCK_METHOD(double, ShapeValue, (const int, const int), (const, override));
  MOCK_METHOD(double, FaceShapeValue, (const int, const int), (const, override));
  MOCK_METHOD((dealii::Tensor<1, dim>), ShapeGradient, (const int, const int), (const, override));
  MOCK_METHOD(double, Jacobian, (const int), (const, override));
  MOCK_METHOD(double, FaceJacobian, (const int), (const, override));
  MOCK_METHOD((dealii::Tensor<1, dim>), FaceNormal, (), (const, override));
  MOCK_METHOD(std::vector<double>, ValueAtQuadrature, (const DealiiVector&), (const, override));
  MOCK_METHOD(std::vector<double>, ValueAtFaceQuadrature, (const DealiiVector&), (const, override));
  MOCK_METHOD(std::vector<Tensor>, GradientAtQuadrature, (const DealiiVector& values_at_dofs), (const, override));
  MOCK_METHOD((dealii::FiniteElement<dim, dim>*), finite_element, (), (override));
  MOCK_METHOD(dealii::FEValues<dim>*, values, (), (override));
  MOCK_METHOD(dealii::FEFaceValues<dim>*, face_values, (), (override));
  MOCK_METHOD(dealii::FEFaceValues<dim>*, neighbor_face_values, (), (override));
};

} // namespace bart::domain::finite_element

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_HPP_
