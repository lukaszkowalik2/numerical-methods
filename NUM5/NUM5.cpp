#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <filesystem>
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

void printFirstAndLastElements(const std::vector<double>& vec, const std::string& name) {
  int n = vec.size();

  std::cout << name << ": [";
  for (int i = 0; i < std::min(5, n); ++i) {
    std::cout << vec[i];
    if (i < std::min(5, n) - 1) std::cout << ", ";
  }
  std::cout << " ... ";
  for (int i = std::max(0, n - 5); i < n; ++i) {
    std::cout << vec[i];
    if (i < n - 1) std::cout << ", ";
  }
  std::cout << "]\n";
}

std::vector<double> createInitialGuess(int N, double value) {
  return std::vector<double>(N, value);
}

std::pair<std::vector<double>, std::vector<double>> jacobiMethod(const std::vector<std::vector<double>>& A,
                                                                 const std::vector<double>& b,
                                                                 const std::vector<double>& exact_solution,
                                                                 const std::vector<double>& initial_guess,
                                                                 int maxIter = 1000, double tol = 1e-10) {
  int N = A.size();
  std::vector<double> x = initial_guess;
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

    if (iter == maxIter - 1) {
    }
  }

  return {x, errors};
}

std::pair<std::vector<double>, std::vector<double>> gaussSeidelMethod(const std::vector<std::vector<double>>& A,
                                                                      const std::vector<double>& b,
                                                                      const std::vector<double>& exact_solution,
                                                                      const std::vector<double>& initial_guess,
                                                                      int maxIter, double tol) {
  int N = A.size();
  std::vector<double> x = initial_guess;
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

void ensureDirectoryExists(const std::string& dirPath) {
  if (!std::filesystem::exists(dirPath)) {
    std::filesystem::create_directory(dirPath);
  }
}

void saveErrorsToFile(const std::vector<double>& jacobi_errors,
                      const std::vector<double>& gauss_errors,
                      const std::string& filename) {
  std::ofstream file("assets/" + filename);
  file << std::scientific << std::setprecision(10);
  file << "# Iteration Jacobi Gauss-Seidel\n";
  for (size_t i = 0; i < std::max(jacobi_errors.size(), gauss_errors.size()); ++i) {
    file << i << " ";
    if (i < jacobi_errors.size() && jacobi_errors[i] > 0)
      file << std::log10(jacobi_errors[i]) << " ";
    else
      file << "NaN ";
    if (i < gauss_errors.size() && gauss_errors[i] > 0)
      file << std::log10(gauss_errors[i]) << "\n";
    else
      file << "NaN\n";
  }
  file.close();
}

void createGnuplotScript(const std::string& datafile, const std::string& plotFile, double d) {
  std::ofstream script("assets/plot_script.gnu");
  script << "set terminal svg size 800,600\n";
  script << "set output 'assets/" << plotFile << ".svg'\n";
  script << "set title 'Błąd w kolejnych iteracjach (d = " << d << ")'\n";
  script << "set xlabel 'Iteracja'\n";
  script << "set ylabel 'log10(Błąd)'\n";
  script << "set grid\n";
  script << "plot 'assets/" << datafile << "' using 1:2 title 'Jacobi' with lines lw 2, \\\n";
  script << "     'assets/" << datafile << "' using 1:3 title 'Gauss-Seidel' with lines lw 2\n";
  script.close();

  system("gnuplot assets/plot_script.gnu");
}

int main() {
  ensureDirectoryExists("assets");

  int N = 200;
  std::vector<double> d_values = {0.5, 1.201, 1.5, 3.0, 5.0};
  std::vector<double> initial_values = {-1.0, 0.0, 10.0, 100.0};
  int maxIter = 50000;
  double tol = 1e-6;

  std::ofstream summary("assets/convergence_summary.txt");
  summary << "d\tPunkt startowy\tJacobi zbieżność\tJacobi błąd\tGauss-Seidel zbieżność\tGauss-Seidel błąd\n";

  for (double d : d_values) {
    std::cout << "\n=== Testowanie dla d = " << d << " ===\n";
    auto A = createMatrix(N, d);
    auto b = createVectorB(N);
    auto exact_solution = solveWithEigen(A, b);

    printFirstAndLastElements(exact_solution, "Eigen LU");

    for (double init_val : initial_values) {
      auto initial_guess = createInitialGuess(N, init_val);
      std::cout << "\nPunkt startowy: " << init_val << "\n";

      auto start = std::chrono::high_resolution_clock::now();
      auto [jacobiSolution, jacobi_errors] = jacobiMethod(A, b, exact_solution, initial_guess, maxIter, tol);
      auto end = std::chrono::high_resolution_clock::now();
      auto jacobi_time = std::chrono::duration<double>(end - start).count();

      start = std::chrono::high_resolution_clock::now();
      auto [gaussSolution, gauss_errors] = gaussSeidelMethod(A, b, exact_solution, initial_guess, maxIter, tol);
      end = std::chrono::high_resolution_clock::now();
      auto gauss_time = std::chrono::duration<double>(end - start).count();

      std::string filename = "errors_d" + std::to_string(d) + "_start" + std::to_string(init_val) + ".dat";
      saveErrorsToFile(jacobi_errors, gauss_errors, filename);
      createGnuplotScript(filename, "convergence_d" + std::to_string(d) + "_start" + std::to_string(init_val), d);

      summary << d << "\t" << init_val << "\t"
              << (jacobi_errors.back() < tol ? "Tak" : "Nie") << "\t"
              << jacobi_errors.back() << "\t"
              << (gauss_errors.back() < tol ? "Tak" : "Nie") << "\t"
              << gauss_errors.back() << "\n";

      std::cout << "Metoda Jacobiego:\n";
      printFirstAndLastElements(jacobiSolution, "Rozwiązanie Jacobi");
      std::cout << "  Czas wykonania: " << jacobi_time << " s\n";
      std::cout << "  Błąd końcowy: " << jacobi_errors.back() << "\n";
      std::cout << "  Liczba iteracji: " << jacobi_errors.size() << "\n";

      std::cout << "Metoda Gaussa-Seidela:\n";
      printFirstAndLastElements(gaussSolution, "Rozwiązanie gauss Seidel");
      std::cout << "  Czas wykonania: " << gauss_time << " s\n";
      std::cout << "  Błąd końcowy: " << gauss_errors.back() << "\n";
      std::cout << "  Liczba iteracji: " << gauss_errors.size() << "\n";
    }
  }

  summary.close();
  return 0;
}
