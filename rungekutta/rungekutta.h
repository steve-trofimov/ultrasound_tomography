#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <vector>

class emmiter {
private:
    double _x, _y, _nx, _ny;
    double centerOfCircle;
public:
    emmiter(const double& x, const double& y, const double& centerOfCircle) :
            _x(x), _y(y), centerOfCircle(centerOfCircle) {
        double temp_nx = centerOfCircle - _x;
        double temp_ny = centerOfCircle - _y;
        _nx = temp_nx / sqrt(temp_nx * temp_nx + temp_ny * temp_ny);
        _ny = temp_ny / sqrt(temp_nx * temp_nx + temp_ny * temp_ny);
    }
    double x() { return _x; }
    double y() { return _y; }
    double nx() { return _nx; }
    double ny() { return _ny; }
};

int emitterNumberSearch(std::vector<emmiter> emmiters, Eigen::VectorXd coor) {
    double minDist = std::numeric_limits<double>::max();
    int resultNumber;
    for (int i = 0; i < emmiters.size(); i++) {
        double tempDist = std::sqrt((emmiters[i].x() - coor.x()) * (emmiters[i].x() - coor.x()) +
                                    (emmiters[i].y() - coor.y()) * (emmiters[i].y() - coor.y()));
        if (minDist > tempDist) {
            minDist = tempDist;
            resultNumber = i;
        }
    }
    return resultNumber;
}

double distance();

class RungeKutta {
public:
  RungeKutta() {}

  template <class Rhs>
  Eigen::VectorXd oneStep(double t, const Eigen::VectorXd &in, double h,
                          Rhs func) {
    Eigen::VectorXd k1, k2, k3, k4;
    k1 = h * func(t, in);
    k2 = h * func(t + h / 2, in + (h / 2) * k1);
    k3 = h * func(t + h / 2, in + (h / 2) * k2);
    k4 = h * func(t + h, in + h * k3);

    return in + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }

    template <class Rhs>
    std::vector<Eigen::VectorXd> functor(double t0, double tEnd, double tStep,
                 emmiter startEmmiter, Rhs func, double& time) {
        std::vector<Eigen::VectorXd> points;
        Eigen::Vector4d nextStep(startEmmiter.x(), startEmmiter.y(), startEmmiter.nx(), startEmmiter.ny());
        points.push_back(nextStep);

        for (double i = t0; i < tEnd; i += tStep) {

            double distanceToCircle = std::abs(sqrt((nextStep.x() - 1024) * (nextStep.x() - 1024) + (nextStep.y() - 1024) * (nextStep.y() - 1024))
                                               - 1024);
            if (distanceToCircle < 1 && (std::abs(nextStep.x() - startEmmiter.x()) + std::abs(nextStep.y() - startEmmiter.y())) > 5)
                break;
            time = i;

            nextStep = oneStep(i, nextStep, tStep, func);
            points.push_back(nextStep);

        }

        return points;
    }
};