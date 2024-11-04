#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <cmath>
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

void gaussElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& x, std::vector<double>& y, int N) {
  std::vector<std::vector<double>> augmented(N, std::vector<double>(N + 1, 0.0));

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

#endif  // ALGORITHMS_H
