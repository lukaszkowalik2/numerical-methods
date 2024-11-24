#include <Eigen/Dense>
#include <chrono>
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

double calculateError(const std::vector<double>& x1, const std::vector<double>& x2) {
  double error = 0.0;
  for (size_t i = 0; i < x1.size(); ++i) {
    error = std::max(error, std::fabs(x1[i] - x2[i]));
  }
  return error;
}

void printVector(const std::vector<double>& vec, const std::string& name) {
  std::cout << name << ": [";
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i < vec.size() - 1) std::cout << ", ";
  }
  std::cout << "]\n";
}

std::vector<double> jacobiMethod(const std::vector<std::vector<double>>& A, const std::vector<double>& b,
                                 int maxIter = 1000, double tol = 1e-10) {
  int N = A.size();
  std::vector<double> x(N, 0.0);
  std::vector<double> x_new(N);
  double max_diff;

  for (int iter = 0; iter < maxIter; iter++) {
    max_diff = 0.0;

    for (int i = 0; i < N; i++) {
      double sum = b[i];
      for (int j = 0; j < N; j++) {
        if (j != i) {
          sum -= A[i][j] * x[j];
        }
      }
      x_new[i] = sum / A[i][i];
      max_diff = std::max(max_diff, std::abs(x_new[i] - x[i]));
    }

    for (int i = 0; i < N; i++) {
      x[i] = x_new[i];
    }

    if (max_diff <= tol) {
      break;
    }
  }

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
      error = std::max(error, std::fabs(x_new - x[i]));
      x[i] = x_new;
    }
    if (error < tol) break;
  }

  return x;
}

std::vector<double> solveWithEigen(const std::vector<std::vector<double>>& matrix, const std::vector<double>& b) {
  Eigen::MatrixXd A = convertToEigen(matrix);
  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());

  Eigen::VectorXd solution = A.fullPivLu().solve(bEigen);

  return std::vector<double>(solution.data(), solution.data() + solution.size());
}

double measureJacobiTime(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIter, double tol) {
  auto start = std::chrono::high_resolution_clock::now();
  jacobiMethod(A, b, maxIter, tol);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  return elapsed.count();
}

double measureGaussSeidelTime(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIter, double tol) {
  auto start = std::chrono::high_resolution_clock::now();
  gaussSeidelMethod(A, b, maxIter, tol);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  return elapsed.count();
}

int main() {
  int N = 200;
  double d1 = 5.0;
  double d2 = 0.5;
  double d3 = 1.201;

  int maxIter = 50000;
  double tol = 1e-6;

  auto A1 = createMatrix(N, d1);
  auto A2 = createMatrix(N, d2);
  auto A3 = createMatrix(N, d3);
  auto b = createVectorB(N);

  Eigen::VectorXd eigenSolution1 = convertToEigen(A1).fullPivLu().solve(Eigen::VectorXd::Map(b.data(), b.size()));
  std::vector<double> eigenVector1(eigenSolution1.data(), eigenSolution1.data() + eigenSolution1.size());

  double jacobiTime1 = measureJacobiTime(A1, b, maxIter, tol);
  double gaussSeidelTime1 = measureGaussSeidelTime(A1, b, maxIter, tol);

  auto jacobiSolution1 = jacobiMethod(A1, b, maxIter, tol);
  auto gaussSeidelSolution1 = gaussSeidelMethod(A1, b, maxIter, tol);

  std::cout << "==== d = " << d1 << " (silnie dominujące) ====\n";
  std::cout << "Czas Jacobiego: " << jacobiTime1 << " s\n";
  std::cout << "Czas Gaussa-Seidela: " << gaussSeidelTime1 << " s\n";
  std::cout << "Błąd Jacobiego: " << calculateError(jacobiSolution1, eigenVector1) << "\n";
  std::cout << "Błąd Gaussa-Seidela: " << calculateError(gaussSeidelSolution1, eigenVector1) << "\n";

  printVector(jacobiSolution1, "Wynik Jacobiego");
  printVector(gaussSeidelSolution1, "Wynik Gausa");

  auto jacobiSolution2 = jacobiMethod(A2, b, maxIter, tol);
  auto gaussSeidelSolution2 = gaussSeidelMethod(A2, b, maxIter, tol);
  auto eigenSolution2 = solveWithEigen(A2, b);
  std::vector<double> eigenVector2(eigenSolution2.begin(), eigenSolution2.end());

  double jacobiTime2 = measureJacobiTime(A2, b, maxIter, tol);
  double gaussSeidelTime2 = measureGaussSeidelTime(A2, b, maxIter, tol);

  std::cout << "\n==== d = " << d2 << " (brak dominacji) ====\n";
  std::cout << "Czas Jacobiego: " << jacobiTime2 << " s\n";
  std::cout << "Czas Gaussa-Seidela: " << gaussSeidelTime2 << " s\n";
  std::cout << "Błąd Jacobiego: " << calculateError(jacobiSolution2, eigenVector2) << "\n";
  std::cout << "Błąd Gaussa-Seidela: " << calculateError(gaussSeidelSolution2, eigenVector2) << "\n";
  printVector(jacobiSolution2, "Wynik Jacobiego");
  printVector(gaussSeidelSolution2, "Wynik Gausa");

  auto jacobiSolution3 = jacobiMethod(A3, b, maxIter, tol);
  auto gaussSeidelSolution3 = gaussSeidelMethod(A3, b, maxIter, tol);
  auto eigenSolution3 = solveWithEigen(A3, b);
  std::vector<double> eigenVector3(eigenSolution3.begin(), eigenSolution3.end());

  double jacobiTime3 = measureJacobiTime(A3, b, maxIter, tol);
  double gaussSeidelTime3 = measureGaussSeidelTime(A3, b, maxIter, tol);

  std::cout << "\n==== d = " << d3 << " (blisko krytycznego punktu) ====\n";
  std::cout << "Czas Jacobiego: " << jacobiTime3 << " s\n";
  std::cout << "Czas Gaussa-Seidela: " << gaussSeidelTime3 << " s\n";
  std::cout << "Błąd Jacobiego: " << calculateError(jacobiSolution3, eigenVector3) << "\n";
  std::cout << "Błąd Gaussa-Seidela: " << calculateError(gaussSeidelSolution3, eigenVector3) << "\n";
  printVector(jacobiSolution3, "Wynik Jacobiego");
  printVector(gaussSeidelSolution3, "Wynik Gausa");

  return 0;
}
