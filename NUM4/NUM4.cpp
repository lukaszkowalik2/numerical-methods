#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include <numeric>
#include <vector>

std::vector<double> matrixVectorProduct(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vec) {
  std::size_t n = matrix.size();
  std::vector<double> result(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      result[i] += matrix[i][j] * vec[j];
    }
  }
  return result;
}

double sumVector(const std::vector<double> &vec) {
  double sum = 0.0;
  for (const auto &element : vec) {
    sum += element;
  }
  return sum;
}

void printVector(const std::vector<double> &vec) {
  std::cout << "Vector size: " << vec.size() << std::endl;
  for (const auto &element : vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

void printVector2D(const std::vector<std::vector<double>> &vec) {
  std::cout << "Vector2D dimensions: " << vec.size() << " x "
            << (vec.empty() ? 0 : vec[0].size()) << std::endl;

  for (const auto &row : vec) {
    for (const auto &element : row) {
      std::cout << element << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

std::vector<double> solveShermanMorrison(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
  std::size_t size = A.size();
  std::vector<std::vector<double>> matrix(2, std::vector<double>(size));

  for (size_t i = 0; i < size; i++) {
    matrix[0][i] = A[i][i] - 1;
    if (i < size - 1) {
      matrix[1][i] = A[i][i + 1] - 1;
    }
  }

  std::vector<double> q(size, 0.0), w(size, 0.0), y;
  q[size - 1] = b[size - 1] / matrix[0][size - 1];
  w[size - 1] = 1.0 / matrix[0][size - 1];

  for (int i = size - 2; i >= 0; --i) {
    q[i] = (b[i] - matrix[1][i] * q[i + 1]) / matrix[0][i];
    w[i] = (1.0 - matrix[1][i] * w[i + 1]) / matrix[0][i];
  }

  double delta = sumVector(q) / (1.0 + sumVector(w));

  for (std::size_t i = 0; i < size; ++i) {
    y.push_back(q[i] - w[i] * delta);
  }

  return y;
}

std::vector<double> solveWithEigen(size_t N, const std::vector<double> &b) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(N, N, 1.0);
  for (size_t i = 0; i < N; ++i) {
    A(i, i) = 5.0;
    if (i < N - 1) {
      A(i, i + 1) = 3.0;
    }
  }

  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());

  Eigen::VectorXd y = A.fullPivLu().solve(bEigen);

  return std::vector<double>(y.data(), y.data() + y.size());
}

int main() {
  const std::size_t N = 120;

  std::vector<double> b(N, 2.0);

  std::vector<std::vector<double>> A(N, std::vector<double>(N, 1.0));
  for (std::size_t i = 0; i < N; i++) {
    A[i][i] = 5.0;
    if (i < N - 1) {
      A[i][i + 1] = 3.0;
    }
  }

  auto startSherman = std::chrono::high_resolution_clock::now();
  std::vector<double> resultSherman = solveShermanMorrison(A, b);
  auto endSherman = std::chrono::high_resolution_clock::now();

  std::cout << "Rozwiązanie przy użyciu metody Shermana-Morrisona (pierwsze 10 elementów):\n";
  for (std::size_t i = 0; i < 10; ++i) {
    std::cout << resultSherman[i] << " ";
  }
  std::cout << "\nCzas wykonania (Sherman-Morrison): "
            << std::chrono::duration<double>(endSherman - startSherman).count() << " sekund.\n\n";

  auto startEigen = std::chrono::high_resolution_clock::now();
  std::vector<double> resultEigen = solveWithEigen(N, b);
  auto endEigen = std::chrono::high_resolution_clock::now();

  std::cout << "Rozwiązanie przy użyciu Eigen (pierwsze 10 elementów):\n";
  for (std::size_t i = 0; i < 10; ++i) {
    std::cout << resultEigen[i] << " ";
  }
  std::cout << "\nCzas wykonania (Eigen): "
            << std::chrono::duration<double>(endEigen - startEigen).count() << " sekund.\n\n";

  double maxDiff = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    maxDiff = std::max(maxDiff, std::abs(resultSherman[i] - resultEigen[i]));
  }
  std::cout << "Maksymalna różnica między metodami: " << maxDiff << "\n";

  return 0;
}
