#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
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

std::pair<std::vector<double>, std::vector<double>> jacobiMethod(const std::vector<std::vector<double>>& A,
                                                                 const std::vector<double>& b,
                                                                 const std::vector<double>& exact_solution,
                                                                 int maxIter = 1000, double tol = 1e-10) {
  int N = A.size();
  std::vector<double> x(N, 0.0);
  std::vector<double> x_new(N);
  std::vector<double> errors;
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

    errors.push_back(calculateError(x, exact_solution));

    if (max_diff <= tol) {
      break;
    }
  }

  return {x, errors};
}

std::pair<std::vector<double>, std::vector<double>> gaussSeidelMethod(const std::vector<std::vector<double>>& A,
                                                                      const std::vector<double>& b,
                                                                      const std::vector<double>& exact_solution,
                                                                      int maxIter, double tol) {
  int N = A.size();
  std::vector<double> x(N, 0.0);
  std::vector<double> errors;

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

    errors.push_back(calculateError(x, exact_solution));

    if (error < tol) break;
  }

  return {x, errors};
}

std::vector<double> solveWithEigen(const std::vector<std::vector<double>>& matrix, const std::vector<double>& b) {
  Eigen::MatrixXd A = convertToEigen(matrix);
  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());

  Eigen::VectorXd solution = A.fullPivLu().solve(bEigen);

  return std::vector<double>(solution.data(), solution.data() + solution.size());
}

void saveErrorsToFile(const std::vector<double>& jacobi_errors,
                      const std::vector<double>& gauss_errors,
                      const std::string& filename) {
  std::ofstream file(filename);
  file << "# Iteration Jacobi Gauss-Seidel\n";
  for (size_t i = 0; i < std::max(jacobi_errors.size(), gauss_errors.size()); ++i) {
    file << i << " ";
    if (i < jacobi_errors.size())
      file << std::log10(jacobi_errors[i]) << " ";
    else
      file << "- ";
    if (i < gauss_errors.size())
      file << std::log10(gauss_errors[i]) << "\n";
    else
      file << "-\n";
  }
  file.close();
}

void createGnuplotScript(const std::string& datafile, const std::string& plotfile, double d) {
  std::ofstream script("plot_script.gnu");
  script << "set terminal svg\n";
  script << "set output '" << plotfile << ".svg'\n";
  script << "set title 'Błąd w kolejnych iteracjach (d = " << d << ")'\n";
  script << "set xlabel 'Iteracja'\n";
  script << "set ylabel 'log10(Błąd)'\n";
  script << "plot '" << datafile << "' using 1:2 title 'Jacobi' with lines, \\\n";
  script << "     '" << datafile << "' using 1:3 title 'Gauss-Seidel' with lines\n";
  script.close();

  system("gnuplot plot_script.gnu");
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

  // Rozwiązanie dokładne dla wszystkich przypadków
  auto exact_solution1 = solveWithEigen(A1, b);
  auto exact_solution2 = solveWithEigen(A2, b);
  auto exact_solution3 = solveWithEigen(A3, b);

  // Przypadek 1 (d = 5.0)
  auto [jacobiSolution1, jacobi_errors1] = jacobiMethod(A1, b, exact_solution1, maxIter, tol);
  auto [gaussSeidelSolution1, gauss_errors1] = gaussSeidelMethod(A1, b, exact_solution1, maxIter, tol);

  saveErrorsToFile(jacobi_errors1, gauss_errors1, "errors_d1.dat");
  createGnuplotScript("errors_d1.dat", "convergence_d1", d1);

  std::cout << "==== d = " << d1 << " (silnie dominujące) ====\n";
  std::cout << "Błąd końcowy Jacobiego: " << jacobi_errors1.back() << "\n";
  std::cout << "Błąd końcowy Gaussa-Seidela: " << gauss_errors1.back() << "\n";

  // Przypadek 2 (d = 0.5)
  auto [jacobiSolution2, jacobi_errors2] = jacobiMethod(A2, b, exact_solution2, maxIter, tol);
  auto [gaussSeidelSolution2, gauss_errors2] = gaussSeidelMethod(A2, b, exact_solution2, maxIter, tol);

  saveErrorsToFile(jacobi_errors2, gauss_errors2, "errors_d2.dat");
  createGnuplotScript("errors_d2.dat", "convergence_d2", d2);

  std::cout << "\n==== d = " << d2 << " (brak dominacji) ====\n";
  std::cout << "Błąd końcowy Jacobiego: " << jacobi_errors2.back() << "\n";
  std::cout << "Błąd końcowy Gaussa-Seidela: " << gauss_errors2.back() << "\n";

  // Przypadek 3 (d = 1.201)
  auto [jacobiSolution3, jacobi_errors3] = jacobiMethod(A3, b, exact_solution3, maxIter, tol);
  auto [gaussSeidelSolution3, gauss_errors3] = gaussSeidelMethod(A3, b, exact_solution3, maxIter, tol);

  saveErrorsToFile(jacobi_errors3, gauss_errors3, "errors_d3.dat");
  createGnuplotScript("errors_d3.dat", "convergence_d3", d3);

  std::cout << "\n==== d = " << d3 << " (blisko krytycznego punktu) ====\n";
  std::cout << "Błąd końcowy Jacobiego: " << jacobi_errors3.back() << "\n";
  std::cout << "Błąd końcowy Gaussa-Seidela: " << gauss_errors3.back() << "\n";

  return 0;
}
