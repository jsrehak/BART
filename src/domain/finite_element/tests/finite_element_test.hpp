#ifndef BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
#define BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_

#include "domain/finite_element/finite_element.hpp"

#include "test_helpers/test_assertions.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace bart::domain::finite_element::testing {

template <int dim>
class FiniteElementBaseClassTest : public ::testing::Test {
 public:
  FiniteElementBaseClassTest() : dof_handler_(triangulation_) {}
  virtual ~FiniteElementBaseClassTest() = default;
 protected:
  dealii::Triangulation<dim> triangulation_;
  dealii::DoFHandler<dim> dof_handler_;

  void TestSetCell(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestSetCellAndFace(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestValueAtQuadrature(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestGraidentAtQuadrature(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestValueAtFaceQuadrature(domain::finite_element::FiniteElement<dim>* test_fe);
  void SetUp() override {
    dealii::GridGenerator::hyper_cube(triangulation_, -1, 1);
    triangulation_.refine_global(2);
  }
};

template <int dim>
void FiniteElementBaseClassTest<dim>::TestSetCell(FiniteElement<dim> *test_fe) {
  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();
  auto cell_id = cell->id();

  EXPECT_NO_THROW(test_fe->SetCell(cell));
  test_fe->values()->reinit(cell);
  EXPECT_FALSE(test_fe->SetCell(cell)); // Shouldn't change anything
  EXPECT_EQ(cell_id, test_fe->values()->get_cell()->id()); // Cell didn't change

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();

  EXPECT_TRUE(test_fe->SetCell(next_cell));
  // Check changed
  EXPECT_NE(cell_id, test_fe->values()->get_cell()->id());
  EXPECT_EQ(next_cell_id, test_fe->values()->get_cell()->id());
}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestSetCellAndFace(FiniteElement<dim> *test_fe) {
  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();
  auto cell_id = cell->id();
  int face = 0;
  int face_index = cell->face_index(face);

  test_fe->SetFace(cell, domain::FaceIndex(face));

  test_fe->face_values()->reinit(cell, face);

  EXPECT_FALSE(test_fe->SetFace(cell, domain::FaceIndex(face)));
  EXPECT_EQ(cell_id, test_fe->face_values()->get_cell()->id());
  EXPECT_EQ(face_index, test_fe->face_values()->get_face_index());

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();
  int next_face = face + 1;
  int next_face_index = next_cell->face_index(next_face);

  EXPECT_TRUE(test_fe->SetFace(next_cell, domain::FaceIndex(next_face)));
  EXPECT_EQ(next_cell_id, test_fe->face_values()->get_cell()->id());
  EXPECT_NE(face_index, test_fe->face_values()->get_face_index());
  EXPECT_EQ(next_face_index, test_fe->face_values()->get_face_index());
}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestValueAtQuadrature(FiniteElement<dim> *test_fe) {

  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();

  EXPECT_NO_THROW(test_fe->SetCell(cell));

  int n_dofs = dof_handler_.n_dofs();

  std::vector<double> moment_values(n_dofs, 0.5);
  dealii::Vector<double> test_moment(moment_values.begin(), moment_values.end());

  std::vector<double> expected_vector(test_fe->n_cell_quad_pts(), 0.5);

  auto result_vector = test_fe->ValueAtQuadrature(test_moment);

  EXPECT_TRUE(bart::test_helpers::AreEqual(expected_vector, result_vector));
}

template <int dim>
<<<<<<< HEAD
void FiniteElementBaseClassTest<dim>::TestValueAtFaceQuadrature(FiniteElement<dim> *test_fe) {
=======
void FiniteElementBaseClassTest<dim>::TestGraidentAtQuadrature(domain::finite_element::FiniteElement<dim> *test_fe) {
  using Tensor = dealii::Tensor<1, dim>;
  using Vector = dealii::Vector<double>;

  dof_handler_.distribute_dofs(*test_fe->finite_element());
  test_fe->SetCell(dof_handler_.begin_active());

  const int n_dofs{ static_cast<int>(dof_handler_.n_dofs()) };
  const int n_cell_quadrature_points{ test_fe->n_cell_quad_pts() };

  const std::vector<double> constant_vector_values(n_dofs, 1.0);
  const Vector test_vector(constant_vector_values.cbegin(), constant_vector_values.cend());

  auto result_tensor_vector = test_fe->GradientAtQuadrature(test_vector);
  ASSERT_EQ(result_tensor_vector.size(), n_cell_quadrature_points);
  for(const auto& tensor : result_tensor_vector) {
    for (int dir = 0; dir < dim; ++dir)
      EXPECT_NEAR(tensor[dir], 0.0, 1e-6);
  }
}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestValueAtFaceQuadrature(
    FiniteElement<dim> *test_fe) {
>>>>>>> Added FiniteElementI::GradientAtQuadrature.

  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();

  EXPECT_NO_THROW(test_fe->SetFace(cell, domain::FaceIndex(0)));

  dealii::Vector<double> values_at_dofs(dof_handler_.n_dofs());
  values_at_dofs = 0.5;
  std::vector<double> expected_vector(test_fe->n_face_quad_pts(), 0.5);

  auto result_vector = test_fe->ValueAtFaceQuadrature(values_at_dofs);

  EXPECT_TRUE(bart::test_helpers::AreEqual(expected_vector, result_vector));
}

} // namespace bart::domain::finite_element::testing

#endif // BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_