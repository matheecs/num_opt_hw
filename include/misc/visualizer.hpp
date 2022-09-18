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
        plt::xlim(-3, 3);
        plt::ylim(-3, 3);
        plt::set_aspect(1);
        if (save_img) {
            plt::save("traj.png", 300);
        } else {
            plt::show();
        }
    }

    static void visualizeObstacles() {
        const int N = 4;
        struct Vertex {
            double x_;
            double y_;

            Vertex(double x, double y) : x_(x), y_(y) {}
        };
        typedef std::vector<Vertex> ConvexObs;
        ConvexObs obs;
        obs.emplace_back(1, 0);
        obs.emplace_back(0, 1);
        obs.emplace_back(-1, 0);
        obs.emplace_back(0, -1);
        obs.emplace_back(1, 0);

        std::vector<double> x(N + 1), y(N + 1);
        for (int i = 0; i <= N; ++i) {
            x[i] = obs[i].x_;
            y[i] = obs[i].y_;
        }

        plt::plot(x, y, "r");
    }
};

#endif