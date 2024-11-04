#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "algorithms.h"

int main() {
  std::ofstream file("times.csv");
  file << "N,Gauss,Thomas,LU,Eigen\n";

  for (int N = 10; N <= 300; N += 10) {
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    constructMatrix(A, N);
    std::vector<double> x(N);
    for (int i = 0; i < N; ++i) x[i] = i + 1;

    std::vector<double> y_gauss(N, 0.0), y_thomas(N, 0.0), y_bandLU(N, 0.0);

    auto start = std::chrono::high_resolution_clock::now();
    gaussElimination(A, x, y_gauss, N);
    auto end = std::chrono::high_resolution_clock::now();
    double time_gauss = std::chrono::duration<double>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    thomasAlgorithm(A, x, y_thomas, N);
    end = std::chrono::high_resolution_clock::now();
    double time_thomas = std::chrono::duration<double>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    bandLU(A, x, y_bandLU, N);
    end = std::chrono::high_resolution_clock::now();
    double time_lu = std::chrono::duration<double>(end - start).count();

    Eigen::MatrixXd eigenA(N, N);
    Eigen::VectorXd eigenX(N), eigenY(N);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) eigenA(i, j) = A[i][j];
      eigenX(i) = x[i];
    }

    start = std::chrono::high_resolution_clock::now();
    eigenY = eigenA.fullPivLu().solve(eigenX);
    end = std::chrono::high_resolution_clock::now();
    double time_eigen = std::chrono::duration<double>(end - start).count();

    file << N << "," << time_gauss << "," << time_thomas << "," << time_lu << "," << time_eigen << "\n";
  }

  file.close();
  std::cout << "Zapisano wyniki do pliku times.csv" << std::endl;
  return 0;
}
