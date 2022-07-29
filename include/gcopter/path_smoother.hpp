#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include <iomanip>
#include "cubic_spline.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cfloat>
#include <cmath>
#include <iostream>
#include <vector>

namespace path_smoother {

class PathSmoother {
 private:
  cubic_spline::CubicSpline cubSpline;

  int pieceN;
  Eigen::Matrix3Xd diskObstacles;
  double penaltyWeight;
  Eigen::Vector2d headP;
  Eigen::Vector2d tailP;
  Eigen::Matrix2Xd points;
  Eigen::Matrix2Xd gradByPoints;

  lbfgs::lbfgs_parameter_t lbfgs_params;

 private:
  static inline double costFunction(void *ptr, const Eigen::VectorXd &x,
                                    Eigen::VectorXd &g) {
    // TODO
    // std::cout << std::setprecision(4) << "x_current:" << x.transpose()
    //           << std::endl;
    const int n = x.size();
    const int num_inner_points = n / 2;
    Eigen::MatrixXd inner_points;
    inner_points.resize(2, num_inner_points);
    inner_points.row(0) = x.head(num_inner_points);
    inner_points.row(1) = x.tail(num_inner_points);

    // std::cout << "inner_points:\n" << inner_points << std::endl;
    auto instance = reinterpret_cast<path_smoother::PathSmoother *>(ptr);

    instance->cubSpline.setInnerPoints(inner_points);

    Eigen::Matrix2Xd g_raw;
    g_raw.resize(2, num_inner_points);
    g_raw.setZero();

    double cost = 0.0;

    // Add Potential
    Eigen::Matrix2Xd gradByPotential;
    gradByPotential.resize(2, num_inner_points);
    gradByPotential.setZero();

    for (int i = 0; i < num_inner_points; ++i) {
      for (int j = 0; j < instance->diskObstacles.cols(); ++j) {
        Eigen::Vector2d diff =
            inner_points.col(i) - instance->diskObstacles.col(j).head(2);
        double distance = diff.norm();
        double delta = instance->diskObstacles(2, j) - distance;

        if (delta > 0.0) {
          cost += delta;
          gradByPotential.col(i) += (-diff / distance);
        }
      }
    }
    // std::cout << "potential(without weight):" << cost << std::endl;
    // std::cout << "potentialGrad:(wo weight)\n" << gradByPotential <<
    // std::endl;
    cost *= instance->penaltyWeight;

    g_raw += instance->penaltyWeight * gradByPotential;

    // Add Energy
    double energy = 0.0;
    instance->cubSpline.getStretchEnergy(energy);
    // std::cout << "energy:" << energy << std::endl;
    cost += energy;

    Eigen::Matrix2Xd gradByEnergy;
    gradByEnergy.resize(2, instance->pieceN - 1);
    instance->cubSpline.getGrad(gradByEnergy);
    g_raw += gradByEnergy;

    // std::cout << "energyGrad:\n" << gradByEnergy << std::endl;

    // std::cout << "total cost(potential+energy):" << cost << std::endl;
    // std::cout << "total grad\n" << g_raw << std::endl;

    g.head(num_inner_points) = g_raw.row(0).transpose();
    g.tail(num_inner_points) = g_raw.row(1).transpose();

    return cost;
  }

  static int monitorProgress(void *instance, const Eigen::VectorXd &x,
                             const Eigen::VectorXd &g, const double fx,
                             const double step, const int k, const int ls) {
    std::cout << std::setprecision(4)
              << "================================" << std::endl
              << "Iteration: " << k << std::endl
              << "Function Value: " << fx << std::endl
              << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
              << "Variables: " << std::endl
              << x.transpose() << std::endl;
    return 0;
  }

 public:
  inline bool setup(const Eigen::Vector2d &initialP,
                    const Eigen::Vector2d &terminalP, const int &pieceNum,
                    const Eigen::Matrix3Xd &diskObs, const double penaWeight) {
    pieceN = pieceNum;
    diskObstacles = diskObs;
    penaltyWeight = penaWeight;
    headP = initialP;
    tailP = terminalP;

    cubSpline.setConditions(headP, tailP, pieceN);

    points.resize(2, pieceN - 1);
    gradByPoints.resize(2, pieceN - 1);

    std::cout << "diskObstacles:\n" << diskObstacles << std::endl;

    return true;
  }

  inline double optimize(CubicCurve &curve, const Eigen::Matrix2Xd &iniInPs,
                         const double &relCostTol) {
    // TODO
    double finalCost;
    std::cout << "Initial points = \n" << iniInPs << std::endl;
    /* Set the initial guess */
    std::cout << "Set the initial guess" << std::endl;
    Eigen::VectorXd x(Eigen::VectorXd::Zero(2 * (pieceN - 1)));
    x.head(pieceN - 1) = iniInPs.row(0).transpose();
    x.tail(pieceN - 1) = iniInPs.row(1).transpose();

    /* Set the minimization parameters */
    std::cout << "Set the minimization parameters" << std::endl;
    lbfgs::lbfgs_parameter_t params;
    params.g_epsilon = 1.0e-8;
    params.past = 3;
    params.delta = relCostTol;

    /* Start minimization */
    std::cout << "Start minimization" << std::endl;
    int ret = lbfgs::lbfgs_optimize(x, finalCost, costFunction, monitorProgress,
                                    this, params);

    /* Report the result. */
    std::cout << std::setprecision(4)
              << "================================" << std::endl
              << "L-BFGS Optimization Returned: " << ret << std::endl
              << "Minimized Cost: " << finalCost << std::endl
              << "Optimal Variables: " << std::endl
              << x.transpose() << std::endl;

    cubSpline.getCurve(curve);

    return finalCost;
  }
};

}  // namespace path_smoother

#endif
