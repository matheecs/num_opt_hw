#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

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
            const int n = x.size();
            Eigen::VectorXd x_raw = x;
            const int num_inner_points = n / 2;
            Eigen::Map <Eigen::MatrixXd> inner_points(x_raw.data(), 2, num_inner_points);
            PathSmoother *pM = (PathSmoother *) ptr;

            pM->cubSpline.setInnerPoints(inner_points);

            std::cout << "pieceN = " << pM->pieceN << std::endl;

            Eigen::Matrix2Xd g_raw;
            g_raw.resize(2, num_inner_points);
            g_raw.setZero();

            std::cout << "Start add cost and grad" << std::endl;

            // Add Potential
            double cost = 0.0;
            Eigen::Matrix2Xd gradByPotential;
            gradByPotential.resize(2, num_inner_points);
            gradByPotential.setZero();
            std::cout << "before for loop" << std::endl;
            for (int i = 0; i < num_inner_points; ++i) {
                for (int j = 0; j < pM->diskObstacles.cols(); ++j) {
                    Eigen::Vector2d diff = inner_points.col(i) - pM->diskObstacles.col(j).head(2);
                    double distance = diff.norm();
                    double delta = pM->diskObstacles(2, j) - distance;

                    std::cout << "running" << std::endl;

                    if (delta > 0.0) {
                        cost += delta;
                        gradByPotential.col(i) += (-diff / distance);
                    }
                }
            }
            cost *= pM->penaltyWeight;
            g_raw += pM->penaltyWeight * gradByPotential;

            std::cout << "add E" << std::endl;

            // Add Energy
            double energy = 0.0;
            pM->cubSpline.getStretchEnergy(energy);
            cost += energy;

            std::cout << "todo add grad of E" << std::endl;

            Eigen::Matrix2Xd gradByEnergy;
            gradByEnergy.resize(2, pM->pieceN - 1);
            pM->cubSpline.getGrad(gradByEnergy);
            g_raw += gradByEnergy;

            g.head(num_inner_points) = g_raw.row(0).transpose();
            g.tail(num_inner_points) = g_raw.row(1).transpose();

            std::cout << "ok" << std::endl;

            return cost;
        }

        static int monitorProgress(void *instance,
                                   const Eigen::VectorXd &x,
                                   const Eigen::VectorXd &g,
                                   const double fx,
                                   const double step,
                                   const int k,
                                   const int ls) {
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

            return true;
        }

        inline double optimize(CubicCurve &curve, const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol) {
            // TODO
            double finalCost;
            std::cout << iniInPs << std::endl;
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
            int ret = lbfgs::lbfgs_optimize(x,
                                            finalCost,
                                            &costFunction,
                                            monitorProgress,
                                            this,
                                            params);

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
