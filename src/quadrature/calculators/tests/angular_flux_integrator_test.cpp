#include "quadrature/calculators/angular_flux_integrator.hpp"

#include <memory>

#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "quadrature/calculators/tests/angular_flux_integrator_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::ContainerEq;

template <typename DimensionWrapper>
class AngularFluxIntegratorTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using QuadraturePointType = NiceMock<typename quadrature::QuadraturePointMock<dim>>;
  using QuadratureSetType = NiceMock<typename quadrature::QuadratureSetMock<dim>>;
  using TestIntegrator = typename quadrature::calculators::AngularFluxIntegrator<dim>;
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<dealii::Vector<double>>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, VectorPtr>;

  // Test parameters
  static constexpr int n_quadrature_points{ 3 };
  static constexpr int n_total_dofs{ 2 };

  // Test object
  std::unique_ptr<TestIntegrator> test_integrator_{ nullptr };

  // Mock objects
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_{ nullptr };
  std::array<std::shared_ptr<QuadraturePointType>, n_quadrature_points> mock_quadrature_points_;

  // Supporting objects
  VectorMap angular_flux_map_{};
  std::array<Vector, n_total_dofs> expected_net_current_at_dofs;
  std::vector<double> expected_directional_current_at_dofs_{1802.25 * dim, 3604.5 * dim};
  std::vector<double> expected_directional_flux_at_dofs_{ 601.5, 1203 };

  auto SetUp() -> void override;
  auto SetUpExpectedValues() -> void;
};

template <typename DimensionWrapper>
auto AngularFluxIntegratorTest<DimensionWrapper>::SetUp() -> void {
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();
  test_integrator_ = std::move(std::make_unique<TestIntegrator>(quadrature_set_ptr_));

  using Index = quadrature::QuadraturePointIndex;
  const std::array<double, n_quadrature_points> weights{0.15, 0.25, 0.6};
  for (int i = 0; i < n_quadrature_points; ++i) {
    auto quadrature_point_ptr = std::make_shared<QuadraturePointType>();
    mock_quadrature_points_.at(i) = quadrature_point_ptr;
    dealii::Tensor<1, dim> position_tensor;
    double omega_value{ 1 };
    if (i == 1) {
      omega_value = -1;
    } else if (i == 2) {
      omega_value = 2;
    }
    for (int j = 0; j < dim; ++j)
      position_tensor[j] = omega_value;
    ON_CALL(*quadrature_point_ptr, weight()).WillByDefault(Return(weights.at(i)));
    ON_CALL(*quadrature_point_ptr, cartesian_position_tensor()).WillByDefault(Return(position_tensor));
    ON_CALL(*quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillByDefault(Return(quadrature_point_ptr));
    auto vector_ptr = std::make_shared<Vector>(n_total_dofs);
    for (int j = 0; j < n_total_dofs; ++j)
      (*vector_ptr)[j] = std::pow(10.0, i + 1) * (j + 1);
    angular_flux_map_.insert({Index(i), vector_ptr});
  }
  ON_CALL(*quadrature_set_ptr_, size()).WillByDefault(Return(n_quadrature_points));
  SetUpExpectedValues();
}

template <typename DimensionWrapper>
auto AngularFluxIntegratorTest<DimensionWrapper>::SetUpExpectedValues() -> void {
  for (int i = 0; i < n_total_dofs; ++i) {
    Vector expected_vector(dim);
    for (int j = 0; j < dim; ++j)
      expected_vector[j] = 1176.5 * (i + 1);
    expected_net_current_at_dofs.at(i) = expected_vector;
  }
}

TYPED_TEST_SUITE(AngularFluxIntegratorTest, bart::testing::AllDimensions);

TYPED_TEST(AngularFluxIntegratorTest, Constructor) {
  constexpr int dim = this->dim;
  using TestIntegrator = typename quadrature::calculators::AngularFluxIntegrator<dim>;
  using QuadratureSetType = typename quadrature::QuadratureSetMock<dim>;

  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  TestIntegrator integrator(quadrature_set_ptr);
  ASSERT_NE(nullptr, integrator.quadrature_set_ptr());
  EXPECT_EQ(integrator.quadrature_set_ptr(), quadrature_set_ptr.get());
}

TYPED_TEST(AngularFluxIntegratorTest, ConstructorBadDependencies) {
  constexpr int dim = this->dim;
  using TestIntegrator = typename quadrature::calculators::AngularFluxIntegrator<dim>;

  EXPECT_ANY_THROW({
                     TestIntegrator integrator(nullptr);
  });
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalCurrentBadNormal) {
  dealii::Vector<double> normal_vector(this->dim + 1);
  for (int i = 0; i < this->dim + 1; ++i)
    normal_vector[i] = 1.5;
  for (int dof = 0; dof < this->n_total_dofs; ++dof) {
    EXPECT_CALL(*this->quadrature_set_ptr_, size())
        .Times(::testing::AtLeast(0))
        .WillRepeatedly(DoDefault());
    for (int i = 0; i < this->n_quadrature_points; ++i) {
      using Index = quadrature::QuadraturePointIndex;
      EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i)))
      .Times(::testing::AtLeast(0))
      .WillRepeatedly(DoDefault());
    }
    for (auto &quadrature_point : this->mock_quadrature_points_) {
      EXPECT_CALL(*quadrature_point, weight())
          .Times(::testing::AtLeast(0))
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*quadrature_point, cartesian_position_tensor())
          .Times(::testing::AtLeast(0))
          .WillRepeatedly(DoDefault());
    }

    using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
    EXPECT_ANY_THROW({
    [[maybe_unused]] auto result = this->test_integrator_->DirectionalCurrent(
        this->angular_flux_map_, normal_vector, DegreeOfFreedom(dof));
                     });
  }
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalCurrentDoFSpecified) {
  dealii::Vector<double> normal_vector(this->dim);
  for (int i = 0; i < this->dim; ++i)
    normal_vector[i] = 1.5;
  for (int dof = 0; dof < this->n_total_dofs; ++dof) {
    EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());
    for (int i = 0; i < this->n_quadrature_points; ++i) {
      using Index = quadrature::QuadraturePointIndex;
      EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillOnce(DoDefault());
    }
    for (auto &quadrature_point : this->mock_quadrature_points_) {
      EXPECT_CALL(*quadrature_point, weight()).WillOnce(DoDefault());
      EXPECT_CALL(*quadrature_point, cartesian_position_tensor()).WillOnce(DoDefault());
    }

    using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
    auto result = this->test_integrator_->DirectionalCurrent(
        this->angular_flux_map_, normal_vector, DegreeOfFreedom(dof));
    EXPECT_DOUBLE_EQ(result, this->expected_directional_current_at_dofs_.at(dof));
  }
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalCurrent) {
  dealii::Vector<double> normal_vector(this->dim);
  for (int i = 0; i < this->dim; ++i)
    normal_vector[i] = 1.5;

  EXPECT_CALL(*this->quadrature_set_ptr_, size())
      .Times(this->n_total_dofs)
      .WillRepeatedly(DoDefault());
  for (int i = 0; i < this->n_quadrature_points; ++i) {
    using Index = quadrature::QuadraturePointIndex;
    EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i)))
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
  }
  for (auto &quadrature_point : this->mock_quadrature_points_) {
    EXPECT_CALL(*quadrature_point, weight())
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
    EXPECT_CALL(*quadrature_point, cartesian_position_tensor())
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
  }

  using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
  auto result = this->test_integrator_->DirectionalCurrent(this->angular_flux_map_, normal_vector);
  EXPECT_TRUE(bart::test_helpers::AreEqual(result, this->expected_directional_current_at_dofs_));
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalFluxBadNormal) {
  dealii::Vector<double> normal_vector(this->dim + 1);
  for (int i = 0; i < this->dim + 1; ++i)
    normal_vector[i] = 1.5;
  for (int dof = 0; dof < this->n_total_dofs; ++dof) {
    EXPECT_CALL(*this->quadrature_set_ptr_, size())
        .Times(::testing::AtLeast(0))
        .WillRepeatedly(DoDefault());
    for (int i = 0; i < this->n_quadrature_points; ++i) {
      using Index = quadrature::QuadraturePointIndex;
      EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i)))
          .Times(::testing::AtLeast(0))
          .WillRepeatedly(DoDefault());
    }
    for (auto &quadrature_point : this->mock_quadrature_points_) {
      EXPECT_CALL(*quadrature_point, weight())
          .Times(::testing::AtLeast(0))
          .WillRepeatedly(DoDefault());
      EXPECT_CALL(*quadrature_point, cartesian_position_tensor())
          .Times(::testing::AtLeast(0))
          .WillRepeatedly(DoDefault());
    }

    using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
    EXPECT_ANY_THROW({
                       [[maybe_unused]] auto result = this->test_integrator_->DirectionalFlux(
                           this->angular_flux_map_, normal_vector, DegreeOfFreedom(dof));
                     });
  }
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalFlux) {
  dealii::Vector<double> normal_vector(this->dim);
  for (int i = 0; i < this->dim; ++i)
    normal_vector[i] = 1.5;
  EXPECT_CALL(*this->quadrature_set_ptr_, size())
      .Times(this->n_total_dofs)
      .WillRepeatedly(DoDefault());
  for (int i = 0; i < this->n_quadrature_points; ++i) {
    using Index = quadrature::QuadraturePointIndex;
    EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i)))
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
  }
  for (auto &quadrature_point : this->mock_quadrature_points_) {
    EXPECT_CALL(*quadrature_point, weight())
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
    EXPECT_CALL(*quadrature_point, cartesian_position_tensor())
        .Times(this->n_total_dofs)
        .WillRepeatedly(DoDefault());
  }

  using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
  auto result = this->test_integrator_->DirectionalFlux(this->angular_flux_map_, normal_vector);
  EXPECT_TRUE(bart::test_helpers::AreEqual(result, this->expected_directional_flux_at_dofs_));
}

TYPED_TEST(AngularFluxIntegratorTest, DirectionalFluxDofSpecified) {
  dealii::Vector<double> normal_vector(this->dim);
  for (int i = 0; i < this->dim; ++i)
    normal_vector[i] = 1.5;
  for (int dof = 0; dof < this->n_total_dofs; ++dof) {
    EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());
    for (int i = 0; i < this->n_quadrature_points; ++i) {
      using Index = quadrature::QuadraturePointIndex;
      EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillOnce(DoDefault());
    }
    for (auto &quadrature_point : this->mock_quadrature_points_) {
      EXPECT_CALL(*quadrature_point, weight()).WillOnce(DoDefault());
      EXPECT_CALL(*quadrature_point, cartesian_position_tensor()).WillOnce(DoDefault());
    }

    using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
    auto result = this->test_integrator_->DirectionalFlux(this->angular_flux_map_,
                                                          normal_vector,
                                                          DegreeOfFreedom(dof));
    EXPECT_DOUBLE_EQ(result, this->expected_directional_flux_at_dofs_.at(dof));
  }
}

TYPED_TEST(AngularFluxIntegratorTest, EddingtonCurrent) {
  using Vector = dealii::Vector<double>;
  using AngleIndex = quadrature::QuadraturePointIndex;
  constexpr int dim = this->dim;
  const double sigma_t{ 4.0 };
  const int angular_flux_max_index{ 2 }; // Just using the non-specific "index" here, this could be quadrature point or dof
  std::map<AngleIndex, std::vector<Vector>> angular_flux_gradient_map;

  EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());
  for (auto &quadrature_point : this->mock_quadrature_points_) {
    EXPECT_CALL(*quadrature_point, weight()).WillOnce(DoDefault());
    EXPECT_CALL(*quadrature_point, cartesian_position_tensor()).WillOnce(DoDefault());
  }

  for (int angle = 0; angle < this->n_quadrature_points; ++angle) {
    EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(AngleIndex(angle))).WillOnce(DoDefault());

    std::vector<Vector> angle_gradient_vector;
    for (int index = 0; index < angular_flux_max_index; ++index) {
      Vector angular_flux_gradient(dim);
      for (int dir = 0; dir < dim; ++dir)
        angular_flux_gradient[dir] = angle + index + 1;
      angle_gradient_vector.push_back(angular_flux_gradient);
    }
    angular_flux_gradient_map[AngleIndex(angle)] = angle_gradient_vector;
  }

  std::vector<Vector> expected_eddington_current;
  for (int index = 0; index < angular_flux_max_index; ++index)
    expected_eddington_current.push_back(Vector(dim));
  const double eddington_current_at_0{ dim/4.0 * (7.85) };
  const double eddington_current_at_1{ dim/4.0 * (10.65) };
  for (int dir = 0; dir < dim; ++dir) {
    expected_eddington_current.at(0)[dir] = eddington_current_at_0;
    expected_eddington_current.at(1)[dir] = eddington_current_at_1;
  }

  const auto returned_eddington_current = this->test_integrator_->EddingtonCurrent(angular_flux_gradient_map, sigma_t);
  ASSERT_EQ(returned_eddington_current.size(), angular_flux_max_index);
  for (int index = 0; index < angular_flux_max_index; ++index) {
    auto& returned_current = returned_eddington_current.at(index);
    ASSERT_EQ(returned_current.size(), dim);
    EXPECT_TRUE(test_helpers::AreEqual(returned_current, expected_eddington_current.at(index)));
  }
}

TYPED_TEST(AngularFluxIntegratorTest, NetCurrentDofSpecified) {
  for (int dof = 0; dof < this->n_total_dofs; ++dof) {
    EXPECT_CALL(*this->quadrature_set_ptr_, size()).WillOnce(DoDefault());
    for (int i = 0; i < this->n_quadrature_points; ++i) {
      using Index = quadrature::QuadraturePointIndex;
      EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i))).WillOnce(DoDefault());
    }
    for (auto &quadrature_point : this->mock_quadrature_points_) {
      EXPECT_CALL(*quadrature_point, weight()).WillOnce(DoDefault());
      EXPECT_CALL(*quadrature_point, cartesian_position_tensor()).WillOnce(DoDefault());
    }

    using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
    auto result = this->test_integrator_->NetCurrent(this->angular_flux_map_, DegreeOfFreedom(dof));
    EXPECT_EQ(result, this->expected_net_current_at_dofs.at(dof));
  }
}

TYPED_TEST(AngularFluxIntegratorTest, NetCurrent) {
  const int n_total_dofs{ this->n_total_dofs };
  EXPECT_CALL(*this->quadrature_set_ptr_, size())
      .Times(n_total_dofs)
      .WillRepeatedly(DoDefault());
  for (int i = 0; i < this->n_quadrature_points; ++i) {
    using Index = quadrature::QuadraturePointIndex;
    EXPECT_CALL(*this->quadrature_set_ptr_, GetQuadraturePoint(Index(i)))
        .Times(n_total_dofs)
        .WillRepeatedly(DoDefault());
  }
  for (auto &quadrature_point : this->mock_quadrature_points_) {
    EXPECT_CALL(*quadrature_point, weight())
        .Times(n_total_dofs)
        .WillRepeatedly(DoDefault());
    EXPECT_CALL(*quadrature_point, cartesian_position_tensor())
        .Times(n_total_dofs)
        .WillRepeatedly(DoDefault());
  }

  using DegreeOfFreedom = typename quadrature::calculators::AngularFluxIntegrator<this->dim>::DegreeOfFreedom;
  auto result = this->test_integrator_->NetCurrent(this->angular_flux_map_);
  EXPECT_EQ(result.size(), this->expected_net_current_at_dofs.size());
  for (std::size_t i = 0; i < result.size(); ++i)
    EXPECT_EQ(result.at(i), this->expected_net_current_at_dofs.at(i));
}







} // namespace
