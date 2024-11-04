#include <Eigen/Dense>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "algorithms.h"

int main() {
  int N = 300;
  std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
  constructMatrix(A, N);

  std::vector<double> x(N);
  for (int i = 0; i < N; ++i) {
    x[i] = i + 1;
  }

  auto myDeterminant = calculateDeterminant(A, N);

  std::vector<double> y_gauss(N, 0.0), y_thomas(N, 0.0), y_bandLU(N, 0.0);

  auto start = std::chrono::high_resolution_clock::now();
  gaussElimination(A, x, y_gauss, N);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Czas Gauss: " << std::chrono::duration<double>(end - start).count() << " sekund\n";

  start = std::chrono::high_resolution_clock::now();
  thomasAlgorithm(A, x, y_thomas, N);
  end = std::chrono::high_resolution_clock::now();
  std::cout << "Czas Thomas: " << std::chrono::duration<double>(end - start).count() << " sekund\n";

  start = std::chrono::high_resolution_clock::now();
  bandLU(A, x, y_bandLU, N);
  end = std::chrono::high_resolution_clock::now();
  std::cout << "Czas LU dla pasmowej: " << std::chrono::duration<double>(end - start).count() << " sekund\n";

  Eigen::MatrixXd eigenA(N, N);
  Eigen::VectorXd eigenX(N), eigenY(N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      eigenA(i, j) = A[i][j];
    }
    eigenX(i) = x[i];
  }

  auto eigenDeterminant = eigenA.determinant();

  start = std::chrono::high_resolution_clock::now();
  eigenY = eigenA.fullPivLu().solve(eigenX);
  end = std::chrono::high_resolution_clock::now();
  std::cout << "Czas Eigen: " << std::chrono::duration<double>(end - start).count() << " sekund\n";

  std::cout << "\nWynik (Gauss):\n";
  for (double val : y_gauss) std::cout << val << " ";
  std::cout << "\n\nWynik (Thomas):\n";
  for (double val : y_thomas) std::cout << val << " ";
  std::cout << "\n\nWynik (LU dla pasmowej):\n";
  for (double val : y_bandLU) std::cout << val << " ";
  std::cout << "\n\nWynik (Eigen):\n"
            << eigenY.transpose() << std::endl;

  std::cout << "Mój wyznacznik: " << myDeterminant << ", wyznacznik Eigen: " << eigenDeterminant << std::endl;

  return 0;
}