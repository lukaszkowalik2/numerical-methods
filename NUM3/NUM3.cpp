#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

void constructMatrix(std::vector<std::vector<double>>& A, int N) {
  for (int i = 0; i < N; ++i) {
    A[i][i] = 1.01;

    if (i + 1 < N) {
      A[i + 1][i] = 0.3;
    }

    if (i + 2 < N) {
      A[i][i + 1] = 0.2 / (i + 1);
    }

    if (i + 3 < N) {
      A[i][i + 2] = 0.15 / pow(i + 1, 3);
    }
  }
}

void printMatrixEigen(const Eigen::MatrixXd& matrix) {
  std::cout << "Macierz czterodiagonalna A:" << std::endl;
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      std::cout << std::setw(10) << std::setprecision(4) << std::fixed << matrix(i, j);
    }
    std::cout << std::endl;
  }
}

void printMatrix(const std::vector<std::vector<double>>& matrix) {
  std::cout << "Macierz czterodiagonalna A:" << std::endl;
  for (const auto& row : matrix) {
    for (const auto& elem : row) {
      std::cout << std::setw(10) << std::setprecision(4) << std::fixed << elem;
    }
    std::cout << std::endl;
  }
}

void gaussElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y, int N) {
  std::vector<std::vector<double>> augmented(N, std::vector<double>(N + 1, 0.0));

  // Initialize augmented matrix with matrix A and vector x
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      augmented[i][j] = A[i][j];
    }
    augmented[i][N] = x[i];
  }

  // Forward elimination
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      double factor = augmented[j][i] / augmented[i][i];
      for (int k = i; k <= N; ++k) {
        augmented[j][k] -= factor * augmented[i][k];
      }
    }
  }

  // Back substitution
  for (int i = N - 1; i >= 0; --i) {
    y[i] = augmented[i][N];
    for (int j = i + 1; j < N; ++j) {
      y[i] -= augmented[i][j] * y[j];
    }
    y[i] /= augmented[i][i];
  }
}

void thomasAlgorithm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y, int N) {
  std::vector<double> a(N, 0.0);  // Sub-diagonal (below the main diagonal)
  std::vector<double> b(N, 0.0);  // Main diagonal
  std::vector<double> c(N, 0.0);  // First super-diagonal (above the main diagonal)
  std::vector<double> d(N, 0.0);  // Second super-diagonal
  std::vector<double> xCopy(x);   // Right-hand side (copy of vector x)

  for (int i = 0; i < N; ++i) {
    b[i] = A[i][i];                     // Main diagonal
    if (i > 0) a[i] = A[i][i - 1];      // Sub-diagonal
    if (i < N - 1) c[i] = A[i][i + 1];  // First super-diagonal
    if (i < N - 2) d[i] = A[i][i + 2];  // Second super-diagonal
  }

  // Forward elimination
  for (int i = 1; i < N; ++i) {
    double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    xCopy[i] -= m * xCopy[i - 1];

    if (i < N - 1) {
      c[i] -= m * d[i - 1];
    }
  }

  for (int i = 2; i < N; ++i) {
    double m = d[i - 2] / b[i - 2];
    b[i] -= m * c[i - 1];
    xCopy[i] -= m * xCopy[i - 2];
  }

  // Back substitution
  y.resize(N, 0.0);
  y[N - 1] = xCopy[N - 1] / b[N - 1];
  y[N - 2] = (xCopy[N - 2] - c[N - 2] * y[N - 1]) / b[N - 2];
  for (int i = N - 3; i >= 0; --i) {
    y[i] = (xCopy[i] - c[i] * y[i + 1] - d[i] * y[i + 2]) / b[i];
  }
}

void bandLU(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y, int N) {
  std::vector<std::vector<double>> L(N, std::vector<double>(N, 0.0)), U(N, std::vector<double>(N, 0.0));
  for (int i = 0; i < N; ++i) L[i][i] = 1.0;

  for (int i = 0; i < N; ++i) {
    for (int j = i; j < std::min(i + 3, N); ++j) {
      double sum = 0.0;
      for (int k = 0; k < i; ++k) {
        sum += L[i][k] * U[k][j];
      }
      U[i][j] = A[i][j] - sum;
    }

    for (int j = i + 1; j < std::min(i + 3, N); ++j) {
      double sum = 0.0;
      for (int k = 0; k < i; ++k) {
        sum += L[j][k] * U[k][i];
      }
      L[j][i] = (A[j][i] - sum) / U[i][i];
    }
  }

  std::vector<double> z(N, 0.0);
  for (int i = 0; i < N; ++i) {
    z[i] = x[i];
    for (int j = 0; j < i; ++j) {
      z[i] -= L[i][j] * z[j];
    }
  }

  for (int i = N - 1; i >= 0; --i) {
    y[i] = z[i];
    for (int j = i + 1; j < N; ++j) {
      y[i] -= U[i][j] * y[j];
    }
    y[i] /= U[i][i];
  }
}

double calculateDeterminant(const std::vector<std::vector<double>>& A, int n) {
  std::vector<double> a(n, 0.0);  // Sub-diagonal (below the main diagonal)
  std::vector<double> b(n, 0.0);  // Main diagonal
  std::vector<double> c(n, 0.0);  // First super-diagonal (above the main diagonal)
  std::vector<double> d(n, 0.0);  // Second super-diagonal

  for (int i = 0; i < n; ++i) {
    b[i] = A[i][i];                     // Main diagonal
    if (i > 0) a[i] = A[i][i - 1];      // Sub-diagonal
    if (i < n - 1) c[i] = A[i][i + 1];  // First super-diagonal
    if (i < n - 2) d[i] = A[i][i + 2];  // Second super-diagonal
  }

  // LU decomposition
  for (int i = 1; i < n - 1; ++i) {
    a[i] /= b[i - 1];
    b[i] -= a[i] * c[i - 1];
    c[i] -= a[i] * d[i - 1];
  }
  a[n - 1] /= b[n - 2];
  b[n - 1] -= a[n - 1] * c[n - 2];

  double determinant = 1.0;
  for (int i = 0; i < n; ++i) {
    determinant *= b[i];
  }
  return determinant;
}

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