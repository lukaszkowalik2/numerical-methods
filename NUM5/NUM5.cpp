#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

std::vector<std::vector<double>> createMatrix(int N, double d) {
  std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
  for (int i = 0; i < N; ++i) {
    A[i][i] = d;
    if (i > 0) A[i][i - 1] = 0.5;
    if (i < N - 1) A[i][i + 1] = 0.5;
    if (i > 1) A[i][i - 2] = 0.1;
    if (i < N - 2) A[i][i + 2] = 0.1;
  }
  return A;
}

std::vector<double> createVectorB(int N) {
  std::vector<double> b(N);
  for (int i = 0; i < N; ++i) {
    b[i] = i + 1;
  }
  return b;
}

void printVector(const std::vector<double>& vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i < vec.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << std::endl;
}

Eigen::MatrixXd convertToEigen(const std::vector<std::vector<double>>& matrix) {
  int rows = matrix.size();
  int cols = matrix[0].size();
  Eigen::MatrixXd eigenMatrix(rows, cols);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      eigenMatrix(i, j) = matrix[i][j];
    }
  }
  return eigenMatrix;
}

std::vector<double> jacobiMethod(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIter, double tol) {
  int N = A.size();
  std::vector<double> x(N, 0.0), x_new(N, 0.0);

  for (int iter = 0; iter < maxIter; ++iter) {
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int j = 0; j < N; ++j) {
        if (j != i) {
          sum += A[i][j] * x[j];
        }
      }
      x_new[i] = (b[i] - sum) / A[i][i];
    }

    double error = 0.0;
    for (int i = 0; i < N; ++i) {
      error += std::fabs(x_new[i] - x[i]);
    }
    if (error < tol) break;

    x = x_new;
  }
  std::cout << "Wynik z metodą Jacobiego " << std::endl;
  printVector(x);

  return x;
}

std::vector<double> gaussSeidelMethod(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIter, double tol) {
  int N = A.size();
  std::vector<double> x(N, 0.0);

  for (int iter = 0; iter < maxIter; ++iter) {
    double error = 0.0;
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int j = 0; j < N; ++j) {
        if (j != i) {
          sum += A[i][j] * x[j];
        }
      }
      double x_new = (b[i] - sum) / A[i][i];
      error += std::fabs(x_new - x[i]);
      x[i] = x_new;
    }
    if (error < tol) break;
  }
  std::cout << "Wynik z metodą Gauss Seidel: " << std::endl;
  printVector(x);

  return x;
}

std::vector<double> solveWithEigen(const std::vector<std::vector<double>>& matrix, const std::vector<double>& b) {
  const Eigen::MatrixXd A = convertToEigen(matrix);

  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());

  Eigen::VectorXd y = A.fullPivLu().solve(bEigen);

  std::vector<double> x = std::vector<double>(y.data(), y.data() + y.size());

  std::cout << "Wynik z metodą Eigen: " << std::endl;
  printVector(x);

  return x;
}

int main() {
  int N = 200;

  double d1 = 2.0;
  double d2 = 0.8;
  double d3 = 1.1;

  int maxIter = 1000;
  double tol = 1e-6;

  auto A1 = createMatrix(N, d1);
  auto A2 = createMatrix(N, d2);
  auto A3 = createMatrix(N, d3);
  auto b = createVectorB(N);

  // std::vector<double> exactSolution(N);
  // for (int i = 0; i < N; ++i) {
  //   exactSolution[i] = (i + 1) / d;
  // }

  auto jacobiSolution1 = jacobiMethod(A1, b, maxIter, tol);
  auto gaussSeidelSolution1 = gaussSeidelMethod(A1, b, maxIter, tol);
  auto eigenSolution1 = solveWithEigen(A1, b);
  std::cout << "-----------------------------------" << std::endl;
  auto jacobiSolution2 = jacobiMethod(A2, b, maxIter, tol);
  auto gaussSeidelSolution2 = gaussSeidelMethod(A2, b, maxIter, tol);
  auto eigenSolution2 = solveWithEigen(A2, b);
  std::cout << "-----------------------------------" << std::endl;
  auto jacobiSolution3 = jacobiMethod(A3, b, maxIter, tol);
  auto gaussSeidelSolution3 = gaussSeidelMethod(A3, b, maxIter, tol);
  auto eigenSolution3 = solveWithEigen(A3, b);

  // std::cout << "Różnice względem dokładnego rozwiązania:\n";
  // std::cout << "Index\tJacobi\t\tGauss-Seidel\n";
  // for (int i = 0; i < N; ++i) {
  //   std::cout << i << "\t"
  //             << std::fabs(jacobiSolution[i] - exactSolution[i]) << "\t"
  //             << std::fabs(gaussSeidelSolution[i] - exactSolution[i]) << "\n";
  // }

  return 0;
}
