#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <vector>

#include "cubic_curve.hpp"

namespace cubic_spline {
// clang-format off

    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    class BandedSystem {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q) {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int actualSize = N * (lowerBw + upperBw + 1);
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            return;
        }

        inline void destroy() {
            if (ptrData != nullptr) {
                delete[] ptrData;
                ptrData = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void) {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j) {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorizeLU() {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; ++k) {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; ++i) {
                    if (operator()(i, k) != 0.0) {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; ++j) {
                    cVl = operator()(k, j);
                    if (cVl != 0.0) {
                        for (int i = k + 1; i <= iM; ++i) {
                            if (operator()(i, k) != 0.0) {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template<typename EIGENMAT>
        inline void solve(EIGENMAT &b) const {
            int iM;
            for (int j = 0; j <= N - 1; ++j) {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; ++i) {
                    if (operator()(i, j) != 0.0) {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j) {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; ++i) {
                    if (operator()(i, j) != 0.0) {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template<typename EIGENMAT>
        inline void solveAdj(EIGENMAT &b) const {
            int iM;
            for (int j = 0; j <= N - 1; ++j) {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; ++i) {
                    if (operator()(j, i) != 0.0) {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j) {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; ++i) {
                    if (operator()(j, i) != 0.0) {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
        }
    };

// clang-format on

class CubicSpline {
 public:
  CubicSpline() = default;

  ~CubicSpline() { A.destroy(); }

 private:
  int N;
  Eigen::Vector2d headP;
  Eigen::Vector2d tailP;
  BandedSystem A;
  Eigen::MatrixX2d b;
  Eigen::Matrix2Xd inner_points;
  std::vector<Eigen::Matrix<double, 2, 4>> cMats;
  Eigen::MatrixXd A_inv;
  Eigen::MatrixXd M_c_x_first, M_d_x_first, M_D_x;

 public:
  inline void setConditions(const Eigen::Vector2d &headPos,
                            const Eigen::Vector2d &tailPos,
                            const int &pieceNum) {
    // TODO
    headP = headPos;
    tailP = tailPos;
    N = pieceNum;

    std::cout << "N = " << N << std::endl;

    A.create(N - 1, 3, 3);
    b.resize(N - 1, 2);
    cMats.resize(N);

    Eigen::MatrixXd A_m;
    A_m.resize(N - 1, N - 1);
    A_m.setZero();
    A_m(0, 0) = 4;
    A_m(0, 1) = 1;
    for (int i = 1; i < N - 2; ++i) {
      A_m(i, i - 1) = 1;
      A_m(i, i) = 4;
      A_m(i, i + 1) = 1;
    }
    A_m(N - 2, N - 3) = 1;
    A_m(N - 2, N - 2) = 4;
    A_inv = A_m.inverse();

    std::cout << "A = \n" << A_m << std::endl;
    std::cout << "A_inv = \n" << A_inv << std::endl;

    M_c_x_first.resize(N, N - 1);
    M_c_x_first.setZero();
    M_c_x_first(0, 0) = 3;
    M_c_x_first(N - 1, N - 2) = -3;
    for (int i = 1; i <= N - 2; ++i) {
      M_c_x_first(i, i - 1) = -3;
      M_c_x_first(i, i) = 3;
    }

    std::cout << "M_c_x_first = \n" << M_c_x_first << std::endl;

    M_d_x_first.resize(N, N - 1);
    M_d_x_first.setZero();
    M_d_x_first(0, 0) = -2;
    M_d_x_first(N - 1, N - 2) = 2;
    for (int i = 1; i <= N - 2; ++i) {
      M_d_x_first(i, i - 1) = 2;
      M_d_x_first(i, i) = -2;
    }

    std::cout << "M_d_x_first = \n" << M_d_x_first << std::endl;

    M_D_x.resize(N - 1, N - 1);
    M_D_x.setZero();
    M_D_x(0, 1) = 1;
    M_D_x(N - 2, N - 3) = -1;
    for (int i = 1; i <= N - 3; ++i) {
      M_D_x(i, i - 1) = -1;
      M_D_x(i, i + 1) = 1;
    }
    std::cout << "M_D_x(raw) = \n" << M_D_x << std::endl;
    M_D_x = 3 * A_inv * M_D_x;

    std::cout << "M_D_x = \n" << M_D_x << std::endl;

    return;
  }

  inline void setInnerPoints(const Eigen::Ref<const Eigen::Matrix2Xd> &inPs) {
    inner_points = inPs;
    // TODO
    A.reset();
    b.setZero();

    A(0, 0) = 4;
    A(0, 1) = 1;
    for (int i = 1; i < N - 2; ++i) {
      A(i, i - 1) = 1;
      A(i, i) = 4;
      A(i, i + 1) = 1;
    }
    A(N - 2, N - 3) = 1;
    A(N - 2, N - 2) = 4;

    b.row(0) = 3 * (inPs.col(1).transpose() - headP.transpose());
    for (int i = 1; i < N - 2; ++i) {
      b.row(i) =
          3 * (inPs.col(i + 1).transpose() - inPs.col(i - 1).transpose());
    }
    b.row(N - 2) = 3 * (tailP.transpose() - inPs.col(N - 3).transpose());

    A.factorizeLU();
    A.solve(b);

    Eigen::Matrix<double, 2, 4> coeffMat;

    // i = 0
    coeffMat.col(3) = headP;
    coeffMat.col(2) = Eigen::Vector2d::Zero();
    coeffMat.col(1) = 3 * (inner_points.col(0) - headP) - b.row(0).transpose();
    coeffMat.col(0) = 3 * (headP - inner_points.col(0)) + b.row(0).transpose();
    // d-c-b-a
    cMats[0] = coeffMat;

    // i = 1->(N-2)
    for (int i = 1; i <= N - 2; ++i) {
      coeffMat.col(3) = inner_points.col(i - 1);
      coeffMat.col(2) = b.row(i - 1).transpose();
      coeffMat.col(1) = 3 * (inner_points.col(i) - inner_points.col(i - 1)) -
                        2 * b.row(i - 1).transpose() - b.row(i).transpose();
      coeffMat.col(0) = 2 * (inner_points.col(i - 1) - inner_points.col(i)) +
                        1 * b.row(i - 1).transpose() + b.row(i).transpose();
      cMats[i] = coeffMat;
    }

    // i = N-1
    coeffMat.col(3) = inner_points.col(N - 2);
    coeffMat.col(2) = b.row(N - 2).transpose();
    coeffMat.col(1) =
        3 * (tailP - inner_points.col(N - 2)) - 2 * b.row(N - 2).transpose();
    coeffMat.col(0) =
        2 * (inner_points.col(N - 2) - tailP) + 1 * b.row(N - 2).transpose();
    cMats[N - 1] = coeffMat;

    // std::cout << "cMat BEGIN\n" << std::endl;
    // for (const auto &cMat : cMats) {
    //   std::cout << cMat << std::endl;
    //   std::cout << "===" << std::endl;
    // }
    // std::cout << "cMat END\n" << std::endl;

    return;
  }

  inline void getCurve(CubicCurve &curve) const {
    // TODO
    std::vector<double> durs(N, 1.0);
    auto c = CubicCurve(durs, cMats);
    curve = std::move(c);
    // std::cout << "curve = \n" << curve.getPositions() << std::endl;

    return;
  }

  inline void getStretchEnergy(double &energy) const {
    // TODO
    double E{0};
    for (int i = 0; i < N; ++i) {
      // 4*c*c + 12*c*d + 12*d*d
      E += 4 * cMats[i].col(1).dot(cMats[i].col(1)) +
           12 * cMats[i].col(1).dot(cMats[i].col(0)) +
           12 * cMats[i].col(0).dot(cMats[i].col(0));
    }
    energy = E;
    return;
  }

  inline const Eigen::MatrixX2d &getCoeffs(void) const { return b; }

  inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> gradByPoints) const {
    // TODO
    // Energy -> 8cc'+12c'd+12cd'+24dd'
    gradByPoints.setZero();
    // Piece 0
    Eigen::VectorXd c_diff, d_diff;
    c_diff = M_c_x_first.row(0).transpose() - M_D_x.row(0).transpose();
    d_diff = M_d_x_first.row(0).transpose() + M_D_x.row(0).transpose();

    gradByPoints += 8 * cMats[0].col(1) * c_diff.transpose() +
                    12 * cMats[0].col(0) * c_diff.transpose() +
                    12 * cMats[0].col(1) * d_diff.transpose() +
                    24 * cMats[0].col(0) * d_diff.transpose();
    // Piece 1->(N-2)
    for (int i = 1; i <= N - 2; ++i) {
      c_diff = M_c_x_first.row(i).transpose() -
               2 * M_D_x.row(i - 1).transpose() - M_D_x.row(i).transpose();
      d_diff = M_d_x_first.row(i).transpose() +
               1 * M_D_x.row(i - 1).transpose() + M_D_x.row(i).transpose();
      gradByPoints += 8 * cMats[i].col(1) * c_diff.transpose() +
                      12 * cMats[i].col(0) * c_diff.transpose() +
                      12 * cMats[i].col(1) * d_diff.transpose() +
                      24 * cMats[i].col(0) * d_diff.transpose();
    }
    // Piece (N-1)
    c_diff =
        M_c_x_first.row(N - 1).transpose() - 2 * M_D_x.row(N - 2).transpose();
    d_diff =
        M_d_x_first.row(N - 1).transpose() + 1 * M_D_x.row(N - 2).transpose();
    gradByPoints += 8 * cMats[N - 1].col(1) * c_diff.transpose() +
                    12 * cMats[N - 1].col(0) * c_diff.transpose() +
                    12 * cMats[N - 1].col(1) * d_diff.transpose() +
                    24 * cMats[N - 1].col(0) * d_diff.transpose();
    return;
  }
};

// clang-format off
}

#endif
