#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

#include "gcopter/cubic_curve.hpp"
#include "gcopter/geo_utils.hpp"
#include "gcopter/quickhull.hpp"
#include "gcopter/trajectory.hpp"
#include "misc/matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct Visualization {
  static void visualize(const CubicCurve &curve) {
    if (curve.getPieceNum() > 0) {
      double T = 0.01;
      Eigen::Vector2d lastX = curve.getPos(0.0);
      std::vector<double> x, y;
      x.push_back(lastX(0));
      y.push_back(lastX(1));
      for (double t = T; t < curve.getTotalDuration(); t += T) {
        Eigen::Vector2d X = curve.getPos(t);
        x.push_back(X(0));
        y.push_back(X(1));
      }
      plt::plot(x, y, "r-");
    }

    if (curve.getPieceNum() > 0) {
      Eigen::MatrixXd wps = curve.getPositions();
      const int N = wps.cols();
      std::vector<double> x(N), y(N);
      for (int i = 0; i < N; i++) {
        x.at(i) = wps.col(i)(0);
        y.at(i) = wps.col(i)(1);
      }

      plt::plot(x, y, "bo");
    }
  }

  static void show_or_save(bool save_img = true) {
    plt::xlim(-30, 30);
    plt::ylim(-30, 30);
    plt::set_aspect(1);
    if (save_img) {
      plt::save("traj.png", 300);
    } else {
      plt::show();
    }
  }

  static void visualizeDisks(const Eigen::Matrix3Xd &disks) {
    for (int i = 0; i < disks.cols(); ++i)
      drawCircle(disks(0, i), disks(1, i), disks(2, i));
  }

  static void drawCircle(const double origin_x, const double origin_y,
                         const double radius) {
    const int N = 50;
    std::cout << origin_x << " " << origin_y << " " << radius << std::endl;
    std::vector<double> x(N), y(N);
    for (int i = 0; i < N; ++i) {
      double t = 2 * M_PI * i / N;
      x.at(i) = radius * cos(t) + origin_x;
      y.at(i) = radius * sin(t) + origin_y;
    }
    x.push_back(radius + origin_x);
    y.push_back(origin_y);
    plt::plot(x, y, "r");
  }
};
#endif