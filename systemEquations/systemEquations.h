#include <ceres/cubic_interpolation.h>

class systemEquation {
  Eigen::MatrixXd speedMap;
  double c0;

public:
  systemEquation(const Eigen::MatrixXd &speedMap, const double &c0);
  Eigen::Vector4d operator()(double t, const Eigen::Vector4d &variables);
};

class contSystem {
  double c0;

public:
  contSystem(const double &c0);
  Eigen::Vector4d operator()(double t, const Eigen::Vector4d &variables);
};
