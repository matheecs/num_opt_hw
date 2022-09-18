#include <yaml-cpp/yaml.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "gcopter/cubic_curve.hpp"
#include "gcopter/cubic_spline.hpp"
#include "gcopter/path_smoother.hpp"
#include "misc/backward.hpp"
#include "misc/visualizer.hpp"

namespace backward {
    backward::SignalHandling sh;
}

struct Config {
    double penaltyWeight;
    Eigen::Matrix2Xd convexObs;
    double pieceLength;
    double relCostTol;

    Config(const std::string &config_file) {
        YAML::Node config = YAML::LoadFile(config_file);
        std::vector<double> circleObsVec;

        penaltyWeight = config["PenaltyWeight"].as<double>();
        circleObsVec = config["CircleObs"].as<std::vector<double>>();
        pieceLength = config["PieceLength"].as<double>();
        relCostTol = config["RelCostTol"].as<double>();

        convexObs = Eigen::Map<const Eigen::Matrix<double, 2, -1, Eigen::ColMajor>>(
                circleObsVec.data(), 2, circleObsVec.size() / 2);
    }
};

class CurveGen {
private:
    Config config;

    std::vector<Eigen::Vector2d> startGoal;

    CubicCurve curve;

public:
    CurveGen() : config("config/curve_gen.yaml") {}

    inline void vizObs() {
        Visualization::visualizeObstacles();
    }

    inline void plan() {
        if (startGoal.size() == 2) {
            const int N =
                    (startGoal.back() - startGoal.front()).norm() / config.pieceLength;
            Eigen::Matrix2Xd innerPoints(2, N - 1);
            std::vector<double> x(N - 1), y(N - 1);
            for (int i = 0; i < N - 1; ++i) {
                innerPoints.col(i) =
                        (startGoal.back() - startGoal.front()) * (i + 1.0) / N +
                        startGoal.front();
                x[i] = innerPoints.col(i)(0);
                y[i] = innerPoints.col(i)(1);
            }
            plt::plot(x, y, "g");


            path_smoother::PathSmoother pathSmoother;
            pathSmoother.setup(startGoal.front(), startGoal.back(), N,
                               config.convexObs, config.penaltyWeight);
            CubicCurve curve;

            if (std::isinf(
                    pathSmoother.optimize(curve, innerPoints, config.relCostTol))) {
                return;
            }

            if (curve.getPieceNum() > 0) {
                std::cout << "Visualization" << std::endl;
                Visualization::visualize(curve);
            }
        }
    }

    inline void run() {
        startGoal.emplace_back(-2, 0.0);
        startGoal.emplace_back(+2, 0.5);
        plan();
        return;
    }
};

int main(int argc, char **argv) {
    CurveGen curveGen;
    // Generate obstacles
    curveGen.vizObs();
    // Generate path
    curveGen.run();
    // Show path
    Visualization::show_or_save(false);
    return 0;
}
