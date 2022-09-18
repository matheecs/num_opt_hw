#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include <Eigen/Eigen>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cubic_spline.hpp"
#include "lbfgs.hpp"
#include "sdqp.hpp"

namespace path_smoother {

    class PathSmoother {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        Eigen::Matrix2Xd convexObstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;


        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr, const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g) {
            const int n = x.size();
            const int num_inner_points = n / 2;
            Eigen::MatrixXd inner_points;
            inner_points.resize(2, num_inner_points);
            inner_points.row(0) = x.head(num_inner_points);
            inner_points.row(1) = x.tail(num_inner_points);

            // std::cout << "inner_points:\n" << inner_points << std::endl;
            auto instance = reinterpret_cast<path_smoother::PathSmoother *>(ptr);

            instance->cubSpline.setInnerPoints(inner_points);

            Eigen::Matrix2Xd grad_total;
            grad_total.resize(2, num_inner_points);
            grad_total.setZero();

            double cost_total = 0.0;

            // ************************* Calculate Potential *************************
            Eigen::Matrix2Xd grad_potential;
            grad_potential.resize(2, num_inner_points);
            grad_potential.setZero();

            Eigen::Matrix<double, 4, 2> A;
            Eigen::Vector4d b;
            A << 1, 1, -1, -1, 1, -1, -1, 1;
            b << 1, 1, 1, 1;
            const double safe_distance = 0.1;
            for (int i = 0; i < num_inner_points; ++i) {
                Eigen::Vector4d A_x = A * inner_points.col(i);
                if (A_x(0) <= b(0)
                    && A_x(1) <= b(1)
                    && A_x(2) <= b(2)
                    && A_x(3) <= b(3)) {
                    // Inside Obs
                    double min_distance = std::numeric_limits<double>::max();
                    Eigen::Vector2d closest_point;
                    for (int j = 0; j < 4; ++j) {
                        double d_point2line = std::abs(A.row(j) * inner_points.col(i) - b(j))
                                              / A.row(j).norm();
                        if (d_point2line < min_distance) {
                            min_distance = d_point2line;
                            closest_point(0) =
                                    A.row(j)(1) *
                                    (
                                            A.row(j)(1) * inner_points.col(i)(0)
                                            - A.row(j)(0) * inner_points.col(i)(1)
                                    )
                                    + A.row(j)(0) * b(j);
                            closest_point(1) =
                                    A.row(j)(0) *
                                    (
                                            -A.row(j)(1) * inner_points.col(i)(0)
                                            + A.row(j)(0) * inner_points.col(i)(1)
                                    )
                                    + A.row(j)(1) * b(j);
                            closest_point /= A.row(j).norm();
                        }
                    }

                    cost_total += safe_distance + min_distance;
                    Eigen::Vector2d diff = inner_points.col(i) - closest_point;
                    grad_potential.col(i) += diff / diff.norm();
                } else {
                    // Outside Obs
                    /*
                     *     min 0.5 x' Q x + c' x,
                     *     s.t. A x <= b,
                     */
                    Eigen::Matrix2d Q;
                    Eigen::Vector2d c;
                    Eigen::Vector2d closest_point;
                    Q << 2, 0, 0, 2;
                    c = -2 * inner_points.col(i);
                    sdqp::sdqp<2>(Q, c, A, b, closest_point);
                    double min_distance = (closest_point - inner_points.col(i)).norm();
                    // TODO
                    if (min_distance < safe_distance) {
                        cost_total += safe_distance - min_distance;
                        Eigen::Vector2d diff = inner_points.col(i) - closest_point;
                        grad_potential.col(i) += -diff / diff.norm();
                    }
                }
                /*
                for (int j = 0; j < instance->convexObstacles.cols(); ++j) {
                    Eigen::Vector2d diff =
                            inner_points.col(i) - instance->convexObstacles.col(j).head(2);
                    double distance = diff.norm();
                    double delta = instance->convexObstacles(2, j) - distance;

                    if (delta > 0.0) {
                        cost_total += delta;
                        grad_potential.col(i) += (-diff / distance);
                    }
                }
                */
            }

            // std::cout << "potential(without weight):" << cost_total << std::endl;
            std::cout << "potentialGrad:(wo weight)\n" << grad_potential << std::endl;
            cost_total *= instance->penaltyWeight;
            grad_total += instance->penaltyWeight * grad_potential;

            // ************************* Calculate Energy *************************
            double cost_energy = 0.0;
            instance->cubSpline.getStretchEnergy(cost_energy);
            // std::cout << "cost_energy:" << cost_energy << std::endl;
            cost_total += cost_energy;

            Eigen::Matrix2Xd grad_energy;
            grad_energy.resize(2, instance->pieceN - 1);
            instance->cubSpline.getGrad(grad_energy);
            grad_total += grad_energy;

            std::cout << "grad_energy:\n" << grad_energy << std::endl;
            // std::cout << "total cost_total:" << cost_total << std::endl;
            std::cout << "total grad\n" << grad_total << std::endl;

            g.head(num_inner_points) = grad_total.row(0).transpose();
            g.tail(num_inner_points) = grad_total.row(1).transpose();

            return cost_total;
        }

        static

        int monitorProgress(void *instance, const Eigen::VectorXd &x,
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
                          const Eigen::Matrix2Xd &convexObs, const double penaWeight) {
            pieceN = pieceNum;
            convexObstacles = convexObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;
            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);

            std::cout << "convexObstacles:\n" << convexObstacles << std::endl;
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
