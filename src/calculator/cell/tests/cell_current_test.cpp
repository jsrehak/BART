#include "calculator/cell/cell_current.hpp"

#include "domain/finite_element/tests/finite_element_mock.hpp"
#include "quadrature/calculators/tests/angular_flux_integrator_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"
#include "test_helpers/dealii_test_domain.h"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return;

template <typename DimensionWrapper>
 class CalculatorCellCurrent : public ::testing::Test, public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using CurrentCalculator = calculator::cell::CellCurrent<dim>;
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorMock;
  using FiniteElement = domain::finite_element::FiniteElementMock<dim>;
  using Vector = dealii::Vector<double>;
  using QuadraturePointIndex = quadrature::QuadraturePointIndex;
  using AngularFluxMap = std::map<QuadraturePointIndex, std::shared_ptr<Vector>>;

  // Test object
  std::unique_ptr<CurrentCalculator> test_calculator_;

  // Dependencies
  std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr_{ std::make_shared<AngularFluxIntegrator>() };
  std::shared_ptr<FiniteElement> finite_element_ptr_{ std::make_shared<FiniteElement>() };

  // Test parameters & supporting objects
  AngularFluxMap angular_flux_map_;
  const int n_quadrature_points{ test_helpers::RandomInt(5, 10) };
  const int n_cell_quadrature_points{ test_helpers::RandomInt(3, 6) };
  const int angular_flux_size{ test_helpers::RandomInt(10, 20) };

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto CalculatorCellCurrent<DimensionWrapper>::SetUp() -> void {
  test_calculator_ = std::make_unique<CurrentCalculator>(finite_element_ptr_, angular_flux_integrator_ptr_);

  for (int i = 0; i < n_quadrature_points; ++i)
    angular_flux_map_[QuadraturePointIndex(i)] = std::make_shared<dealii::Vector<double>>(angular_flux_size);

  this->SetUpDealii();
}

TYPED_TEST_SUITE(CalculatorCellCurrent, bart::testing::AllDimensions);

TYPED_TEST(CalculatorCellCurrent, Constructor) {
  constexpr int dim = this->dim;
  using CurrentCalculator = calculator::cell::CellCurrent<dim>;
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorMock;
  using FiniteElement = domain::finite_element::FiniteElementMock<dim>;

  EXPECT_NO_THROW({
    CurrentCalculator test_calculator(std::make_shared<FiniteElement>(), std::make_shared<AngularFluxIntegrator>());
  });
}

TYPED_TEST(CalculatorCellCurrent, ConstructorBadDependency) {
  constexpr int dim = this->dim;
  using CurrentCalculator = calculator::cell::CellCurrent<dim>;
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorMock;
  using FiniteElement = domain::finite_element::FiniteElementMock<dim>;

  EXPECT_ANY_THROW({ CurrentCalculator test_calculator(nullptr, std::make_shared<AngularFluxIntegrator>()); });
  EXPECT_ANY_THROW({ CurrentCalculator test_calculator(std::make_shared<FiniteElement>(), nullptr); });
}

TYPED_TEST(CalculatorCellCurrent, Getters) {
  EXPECT_EQ(this->test_calculator_->finite_element_ptr(), this->finite_element_ptr_.get());
  EXPECT_EQ(this->test_calculator_->angular_flux_integrator_ptr(), this->angular_flux_integrator_ptr_.get());
}

TYPED_TEST(CalculatorCellCurrent, CurrentAtQuadrature) {
  constexpr int dim{ this->dim };
  using Vector = dealii::Vector<double>;
  auto& cell = *this->cells_.begin();

  EXPECT_CALL(*this->finite_element_ptr_, SetCell(cell)).WillOnce(Return(true));

  /* We need to create a vector of vectors, that is the current at each DOF, as calculated by the
   * net current function, based on the angular fluxes given to it. We will also populate Vectors that
   * represent the directional components for the current. Each of these will be passed to the ValueAtQuadrature
   * and then those will be reassembled to get the current at the quadrature points
   *
   * To summarize, currents at dofs is a vector of size dofs of vectors of size dim, and current_components_at_dofs
   * is an array of size dim of vectors of size dofs! */
  std::vector<Vector> currents_at_dofs(this->angular_flux_size);
  std::array<Vector, dim> current_components_at_dofs;
  for (auto& current : currents_at_dofs) {
    const auto net_current_at_dof{ test_helpers::RandomVector(dim, -100, 100) };
    current = Vector(net_current_at_dof.cbegin(), net_current_at_dof.cend());
  }
  for (int dir = 0; dir < dim; ++dir) {
    current_components_at_dofs.at(dir).reinit(this->angular_flux_size);
    for (int i = 0; i < this->angular_flux_size; ++i) {
      current_components_at_dofs.at(dir)[i] = currents_at_dofs.at(i)[dir];
    }
  }

  // This is the current copmonents at the quadrature points. We are doing the opposite we did above
  std::vector<Vector> current_at_quadrature(this->n_cell_quadrature_points);
  std::array<std::vector<double>, dim> current_components_at_quadrature;
  for (int dir = 0; dir < dim; ++dir) {
    auto& current_componant = current_components_at_quadrature.at(dir);
    current_componant.resize(this->n_cell_quadrature_points);
    for (int i = 0; i < this->n_cell_quadrature_points; ++i) {
      current_componant[i] = test_helpers::RandomDouble(-100, 100);
    }
  }
  for(int i = 0; i < this->n_cell_quadrature_points; ++i) {
    current_at_quadrature.at(i).reinit(dim);
    for (int dir = 0; dir < dim; ++dir) {
      current_at_quadrature.at(i)[dir] = current_components_at_quadrature.at(dir).at(i);
    }
  }

  EXPECT_CALL(*this->angular_flux_integrator_ptr_, NetCurrent(this->angular_flux_map_))
      .WillOnce(Return(currents_at_dofs));
  for (int i = 0; i < dim; ++i) {
    EXPECT_CALL(*this->finite_element_ptr_, ValueAtQuadrature(current_components_at_dofs.at(i)))
        .WillOnce(Return(current_components_at_quadrature.at(i)));
  }

  auto returned_current = this->test_calculator_->CurrentAtQuadrature(cell, this->angular_flux_map_);
  ASSERT_EQ(returned_current.size(), this->n_cell_quadrature_points);
  for (int i = 0; i < this->n_cell_quadrature_points; ++i) {
    EXPECT_TRUE(test_helpers::AreEqual(returned_current.at(i), current_at_quadrature.at(i)));
  }
}

} // namespace
