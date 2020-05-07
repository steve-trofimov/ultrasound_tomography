#include "systemEquations.h"

systemEquation::systemEquation(const Eigen::MatrixXd &speedMap,
                               const double &c0)
    : speedMap(speedMap), c0(c0) {}

Eigen::Vector4d systemEquation::operator()(double t,
                                           const Eigen::Vector4d &variables) {
  ceres::Grid2D<double, 1> grid(speedMap.data(), 0, speedMap.rows(), 0,
                                speedMap.cols());
  ceres::BiCubicInterpolator<ceres::Grid2D<double, 1>> interpolator(grid);

  double c, dcdx, dcdy;
  interpolator.Evaluate(variables.x(), variables.y(), &c, &dcdx, &dcdy);

  Eigen::Vector2d n(variables.z(), variables.w());
  n.normalize();

  Eigen::Vector4d result;

  result << c * n.x(), c * n.y(), -c0 * dcdx / c, -c0 * dcdy / c;

  return result;
}

contSystem::contSystem(const double &c0) : c0(c0) {}

Eigen::Vector4d contSystem::operator()(double t,
                                       const Eigen::Vector4d &variables) {
  double c, dcdx, dcdy;
  c = 300. + 500. * exp(-(4 * variables.x() * variables.x() +
                          variables.y() * variables.y()) /
                        (251 * 251));
  dcdx =
      -4000. * variables.x() *
      exp(-(4 * variables.x() * variables.x() + variables.y() * variables.y()) /
          (251 * 251)) /
      251 * 251;
  dcdy =
      -1000. * variables.y() *
      exp(-(4 * variables.x() * variables.x() + variables.y() * variables.y()) /
          (251 * 251)) /
      251 * 251;

  Eigen::Vector2d n(variables.z(), variables.w());
  n.normalize();

  Eigen::Vector4d result;

  result << c * n.x(), c * n.y(), -c0 * dcdx / c, -c0 * dcdy / c;

  return result;
}
