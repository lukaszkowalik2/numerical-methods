#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

void writeConvergenceData(const std::vector<double>& data, const std::string& fname) {
  std::ofstream out(fname + ".dat");
  for (size_t i = 0; i < data.size(); ++i) {
    out << i << " " << std::log10(data[i]) << "\n";
  }
  out.close();
  system("gnuplot plot_convergence.plt");
}

void displayMatrix(const std::vector<std::vector<double>>& matrix, const std::string& header = "") {
  if (!header.empty()) std::cout << header << "\n";

  const int width = 10;
  for (const auto& row : matrix) {
    for (const auto& elem : row) {
      std::cout << std::setw(width) << elem << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void displayVector(const std::vector<double>& vec, const std::string& header = "") {
  if (!header.empty()) std::cout << header << "\n";

  for (const auto& x : vec) {
    std::cout << std::setw(10) << x << " ";
  }
  std::cout << "\n\n";
}

std::vector<double> multiply(const std::vector<std::vector<double>>& mat,
                             const std::vector<double>& vec) {
  std::vector<double> result(4);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      result[i] += mat[i][j] * vec[j];
    }
  }
  return result;
}

std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>>& A,
                                          const std::vector<std::vector<double>>& B) {
  std::vector<std::vector<double>> result(4, std::vector<double>(4));
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return result;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double sum = 0;
  for (int i = 0; i < 4; ++i) sum += a[i] * b[i];
  return sum;
}

double magnitude(const std::vector<double>& v) {
  return std::sqrt(dot(v, v));
}

void powerMethod(const std::vector<std::vector<double>>& matrix,
                 double& eigenvalue,
                 std::vector<double>& eigenvector,
                 std::vector<double>& errors,
                 double tolerance = 1e-12) {
  std::vector<double> x(4, 1.0);
  double norm = magnitude(x);
  for (auto& xi : x) xi /= norm;

  double lambda_prev = 0;
  double lambda = 0;

  while (true) {
    auto y = multiply(matrix, x);
    lambda = dot(x, y);
    errors.push_back(std::abs(lambda - lambda_prev));

    double ny = magnitude(y);
    for (int i = 0; i < 4; ++i) x[i] = y[i] / ny;

    if (std::abs(lambda - lambda_prev) < tolerance) break;
    lambda_prev = lambda;
  }

  eigenvalue = lambda;
  eigenvector = x;
}

void QRDecomposition(const std::vector<std::vector<double>>& A,
                     std::vector<std::vector<double>>& Q,
                     std::vector<std::vector<double>>& R) {
  Q = A;
  for (int j = 0; j < 4; ++j) {
    double col_norm = 0;
    for (int i = 0; i < 4; ++i) col_norm += Q[i][j] * Q[i][j];
    R[j][j] = std::sqrt(col_norm);

    for (int i = 0; i < 4; ++i) Q[i][j] /= (R[j][j] + 1e-30);

    for (int k = j + 1; k < 4; ++k) {
      R[j][k] = 0;
      for (int i = 0; i < 4; ++i) R[j][k] += Q[i][j] * Q[i][k];
      for (int i = 0; i < 4; ++i) Q[i][k] -= R[j][k] * Q[i][j];
    }
  }

  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < i; ++j)
      R[i][j] = 0;
}

double computeWilkinsonShift(const std::vector<std::vector<double>>& A) {
  const int n = A.size();
  const double a = A[n - 2][n - 2];
  const double b = A[n - 2][n - 1];
  const double c = A[n - 1][n - 1];
  const double d = (a - c) / 2.0;
  const double sgn = (d >= 0) ? 1.0 : -1.0;
  return c - (sgn * b * b) / (std::abs(d) + std::sqrt(d * d + b * b));
}

std::vector<std::vector<double>> QRIterationWithShift(std::vector<std::vector<double>> A,
                                                      std::vector<double>& errors,
                                                      double tolerance = 1e-12) {
  double error;
  do {
    const double shift = computeWilkinsonShift(A);
    for (int i = 0; i < 4; ++i) A[i][i] -= shift;

    std::vector<std::vector<double>> Q(4, std::vector<double>(4));
    std::vector<std::vector<double>> R(4, std::vector<double>(4));
    QRDecomposition(A, Q, R);

    A = multiply(R, Q);
    for (int i = 0; i < 4; ++i) A[i][i] += shift;

    error = 0;
    for (int i = 1; i < 4; ++i) error += std::abs(A[i][i - 1]);
    errors.push_back(error);
  } while (error > tolerance);

  return A;
}

void saveComparisonData(const std::vector<double>& basic,
                        const std::vector<double>& wilkinson,
                        const std::string& fname) {
  std::ofstream data(fname + ".dat");
  for (size_t i = 0; i < basic.size(); ++i) {
    data << i << " " << std::log10(basic[i]) << " "
         << std::log10(wilkinson[i]) << "\n";
  }
  data.close();
  system("gnuplot plot_qr_comparison.plt");
}

std::vector<std::vector<double>> basicQRIteration(std::vector<std::vector<double>> A,
                                                  std::vector<double>& errors,
                                                  std::vector<std::vector<std::vector<double>>>& all_iterations,
                                                  double tolerance = 1e-12) {
  all_iterations.push_back(A);
  double error;

  do {
    std::vector<std::vector<double>> Q(4, std::vector<double>(4));
    std::vector<std::vector<double>> R(4, std::vector<double>(4));
    QRDecomposition(A, Q, R);
    A = multiply(R, Q);

    all_iterations.push_back(A);

    error = 0;
    for (int i = 1; i < 4; ++i) error += std::abs(A[i][i - 1]);
    errors.push_back(error);
  } while (error > tolerance);

  return A;
}

void saveDiagonalEvolution(const std::vector<std::vector<std::vector<double>>>& iterations,
                           const std::string& fname) {
  Eigen::Matrix4d eigen_matrix;
  eigen_matrix << 9, 2, 0, 0,
      2, 4, 1, 0,
      0, 1, 3, 1,
      0, 0, 1, 2;

  Eigen::EigenSolver<Eigen::Matrix4d> solver(eigen_matrix);
  std::vector<double> exact_eigenvalues(4);
  for (int i = 0; i < 4; ++i) {
    exact_eigenvalues[i] = solver.eigenvalues()(i).real();
  }

  std::sort(exact_eigenvalues.begin(), exact_eigenvalues.end(), std::greater<double>());

  std::ofstream data(fname + ".dat");

  for (size_t iter = 0; iter < iterations.size(); ++iter) {
    data << iter;
    std::vector<double> current_eigenvalues(4);
    for (int i = 0; i < 4; ++i) {
      current_eigenvalues[i] = iterations[iter][i][i];
    }
    std::sort(current_eigenvalues.begin(), current_eigenvalues.end(), std::greater<double>());

    for (int i = 0; i < 4; ++i) {
      double diff = std::abs(current_eigenvalues[i] - exact_eigenvalues[i]);
      data << " " << (diff > 1e-16 ? std::log10(diff) : -16.0);
    }
    data << "\n";
  }
  data.close();
  system("gnuplot diagonal_evolution.plt");
}

int main() {
  const std::vector<std::vector<double>> matrix = {
      {9, 2, 0, 0},
      {2, 4, 1, 0},
      {0, 1, 3, 1},
      {0, 0, 1, 2}};

  double eigenvalue;
  std::vector<double> eigenvector(4);
  std::vector<double> power_errors;
  powerMethod(matrix, eigenvalue, eigenvector, power_errors);

  std::cout << "Najwieksza (co do modulow) wartosc wlasna: " << eigenvalue << "\n";
  std::cout << "Odpowiadajacy wektor wlasny (znormalizowany):\n";
  displayVector(eigenvector);

  writeConvergenceData(power_errors, "power_method_convergence");
  std::cout << "Wykres zbieznosci zapisano do pliku 'power_method_convergence.png'\n\n";

  std::vector<double> basic_errors;
  std::vector<std::vector<std::vector<double>>> iterations;
  auto final_matrix = basicQRIteration(matrix, basic_errors, iterations);
  displayMatrix(final_matrix, "Macierz po osiagnieciu zbieznosci QR:");

  std::cout << "Wartosci wlasne (odczyt z diagonali macierzy trojkatnej):\n";
  for (int i = 0; i < 4; ++i) {
    std::cout << "lambda_" << i << " ~ " << final_matrix[i][i] << "\n";
  }

  std::cout << "\n=== (b2) Algorytm QR z przesunieciem Wilkinsona ===\n";
  std::vector<double> wilkinson_errors;
  auto wilkinson_matrix = QRIterationWithShift(matrix, wilkinson_errors);
  displayMatrix(wilkinson_matrix, "Macierz po osiagnieciu zbieznosci QR z przesunieciem Wilkinsona:");

  std::cout << "Wartosci wlasne (odczyt z diagonali macierzy):\n";
  for (int i = 0; i < 4; ++i) {
    std::cout << "lambda_" << i << " ~ " << wilkinson_matrix[i][i] << "\n";
  }

  saveComparisonData(basic_errors, wilkinson_errors, "qr_comparison");
  std::cout << "Wykres porownania zbieznosci zapisano do pliku 'qr_comparison.png'\n";

  Eigen::Matrix4d eigen_matrix;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      eigen_matrix(i, j) = matrix[i][j];

  Eigen::EigenSolver<Eigen::Matrix4d> solver(eigen_matrix);
  auto values = solver.eigenvalues();
  auto vectors = solver.eigenvectors();

  std::cout << "Wartosci wlasne z Eigen:\n";
  for (int i = 0; i < 4; ++i) {
    std::cout << "lambda_" << i << " = " << values(i).real()
              << (std::abs(values(i).imag()) > 1e-10 ? " + " + std::to_string(values(i).imag()) + "i" : "")
              << "\n";
  }

  std::cout << "\n=== Analiza zbieznosci do macierzy trojkatnej gornej ===\n";
  std::cout << "Pierwsza iteracja (i=0):\n";
  displayMatrix(iterations[0]);

  int centerElement = iterations.size() / 2;
  std::cout << "Srodkowa iteracja (i=" << centerElement << "):\n";
  displayMatrix(iterations[centerElement]);

  int lastElement = iterations.size() - 1;
  std::cout << "Ostatnia iteracja (i=" << lastElement << "):\n";
  displayMatrix(iterations[lastElement]);

  saveDiagonalEvolution(iterations, "diagonal_evolution");

  return 0;
}
