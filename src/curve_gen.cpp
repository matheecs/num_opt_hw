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
  std::string targetTopic;
  double penaltyWeight;
  Eigen::Matrix3Xd circleObs;
  double pieceLength;
  double relCostTol;

  Config(const std::string &config_file) {
    YAML::Node config = YAML::LoadFile(config_file);
    std::vector<double> circleObsVec;

    targetTopic = config["TargetTopic"].as<std::string>();
    penaltyWeight = config["PenaltyWeight"].as<double>();
    circleObsVec = config["CircleObs"].as<std::vector<double>>();
    pieceLength = config["PieceLength"].as<double>();
    relCostTol = config["RelCostTol"].as<double>();

    circleObs = Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>>(
        circleObsVec.data(), 3, circleObsVec.size() / 3);
  }
};

class CurveGen {
 private:
  Config config;

  std::vector<Eigen::Vector2d> startGoal;

  CubicCurve curve;

 public:
  CurveGen() : config("config/curve_gen.yaml") {}
  inline void vizObs() { Visualization::visualizeDisks(config.circleObs); }

  inline void plan() {
    if (startGoal.size() == 2) {
      const int N =
          (startGoal.back() - startGoal.front()).norm() / config.pieceLength;
      Eigen::Matrix2Xd innerPoints(2, N - 1);
      for (int i = 0; i < N - 1; ++i) {
        innerPoints.col(i) =
            (startGoal.back() - startGoal.front()) * (i + 1.0) / N +
            startGoal.front();
      }

      path_smoother::PathSmoother pathSmoother;
      pathSmoother.setup(startGoal.front(), startGoal.back(), N,
                         config.circleObs, config.penaltyWeight);
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
    startGoal.emplace_back(20, -22);
    startGoal.emplace_back(-6, -6);
    plan();
    return;
  }
};

int main(int argc, char **argv) {
  CurveGen curveGen;
  curveGen.vizObs();
  curveGen.run();
  Visualization::show_or_save(false);
  return 0;
}
