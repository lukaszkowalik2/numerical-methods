#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

double sumVector(const std::vector<double>& vec) {
  return std::accumulate(vec.begin(), vec.end(), 0.0);
}

void printVector(const std::vector<double>& vec, const std::string& name = "") {
  std::cout << "Rozmiar wektora: " << vec.size() << std::endl;
  if (name != "") {
    std::cout << "Rozwiązanie dla " << name << std::endl;
  }
  for (const auto& element : vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

void printFirstTenElements(const std::vector<double>& vec, const std::string& name) {
  std::cout << "Pierwsze 10 elementów " << name << ":" << std::endl;
  const size_t elementsToShow = std::min(size_t(10), vec.size());
  for (size_t i = 0; i < elementsToShow; ++i) {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}

void printVector2D(const std::vector<std::vector<double>>& vec) {
  if (vec.empty()) {
    std::cout << "Pusta macierz" << std::endl;
    return;
  }

  std::cout << "Wymiary macierzy: " << vec.size() << " x " << vec[0].size() << std::endl;
  for (const auto& row : vec) {
    for (const auto& element : row) {
      std::cout << element << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

std::vector<double> solveShermanMorrison(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
  const std::size_t size = A.size();
  std::vector<std::vector<double>> matrix(2, std::vector<double>(size));

  for (size_t i = 0; i < size; i++) {
    matrix[0][i] = A[i][i] - 1;
    if (i < size - 1) {
      matrix[1][i] = A[i][i + 1] - 1;
    }
  }

  std::vector<double> q(size), w(size), y(size);
  q[size - 1] = b[size - 1] / matrix[0][size - 1];
  w[size - 1] = 1.0 / matrix[0][size - 1];

  for (int i = size - 2; i >= 0; --i) {
    q[i] = (b[i] - matrix[1][i] * q[i + 1]) / matrix[0][i];
    w[i] = (1.0 - matrix[1][i] * w[i + 1]) / matrix[0][i];
  }

  const double delta = sumVector(q) / (1.0 + sumVector(w));

  for (std::size_t i = 0; i < size; ++i) {
    y[i] = q[i] - w[i] * delta;
  }

  return y;
}

Eigen::MatrixXd createMatrix(size_t N) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(N, N, 1.0);
  for (size_t i = 0; i < N; ++i) {
    A(i, i) = 5.0;
    if (i < N - 1) {
      A(i, i + 1) = 3.0;
    }
  }
  return A;
}

std::vector<double> solveWithEigen(size_t N, const std::vector<double>& b) {
  Eigen::MatrixXd A = createMatrix(N);
  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());
  Eigen::VectorXd y = A.fullPivLu().solve(bEigen);
  return std::vector<double>(y.data(), y.data() + y.size());
}

std::vector<double> solveWithGauss(size_t N, const std::vector<double>& b) {
  Eigen::MatrixXd A = createMatrix(N);
  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());
  Eigen::VectorXd y = A.colPivHouseholderQr().solve(bEigen);
  return std::vector<double>(y.data(), y.data() + y.size());
}

std::vector<double> solveWithQR(size_t N, const std::vector<double>& b) {
  Eigen::MatrixXd A = createMatrix(N);
  Eigen::VectorXd bEigen = Eigen::VectorXd::Map(b.data(), b.size());
  Eigen::VectorXd y = A.householderQr().solve(bEigen);
  return std::vector<double>(y.data(), y.data() + y.size());
}

template <typename Func>
double measureExecutionTime(Func&& func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(end - start).count();
}

void createGnuplotScript(const std::string& filename, const std::string& title, bool isLarge = false) {
  std::ofstream script(filename);
  script << "set terminal svg\n"
         << "set output '" << (isLarge ? "timing_comparison_large.svg" : "timing_comparison.svg") << "'\n"
         << "set title '" << title << "'\n"
         << "set xlabel 'Rozmiar macierzy'\n"
         << "set ylabel 'Czas wykonania [s]'\n"
         << "set grid\n";

  if (!isLarge) {
    script << "plot 'timing_data.txt' using 1:2 with linespoints title 'Sherman-Morrison', "
           << "'timing_data.txt' using 1:3 with linespoints title 'Eigen', "
           << "'timing_data.txt' using 1:4 with linespoints title 'Gauss', "
           << "'timing_data.txt' using 1:5 with linespoints title 'QR'\n";
  } else {
    script << "plot 'timing_data_large.txt' using 1:2 with linespoints title 'Sherman-Morrison'\n";
  }
}

double calculateError(const std::vector<double>& x1, const std::vector<double>& x2) {
  if (x1.size() != x2.size()) {
    throw std::runtime_error("Wektory muszą mieć ten sam rozmiar!");
  }

  double maxError = 0.0;
  for (size_t i = 0; i < x1.size(); ++i) {
    double error = std::abs(x1[i] - x2[i]);
    maxError = std::max(maxError, error);
  }
  return maxError;
}

int main() {
  std::ofstream dataFile("timing_data.txt");
  std::ofstream dataFileLarge("timing_data_large.txt");

  const std::vector<size_t> sizes = {10, 20, 60, 80, 100, 120};
  std::vector<size_t> largerSizes;
  for (size_t N = 500; N <= 30000; N += 500) {
    largerSizes.push_back(N);
  }

  for (const auto& N : sizes) {
    std::vector<double> b(N, 2.0);
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 1.0));

    for (std::size_t i = 0; i < N; i++) {
      A[i][i] = 5.0;
      if (i < N - 1) {
        A[i][i + 1] = 3.0;
      }
    }

    std::vector<double> resultSherman;
    double timeSherman = measureExecutionTime([&]() {
      resultSherman = solveShermanMorrison(A, b);
    });

    std::vector<double> resultEigen;
    double timeEigen = measureExecutionTime([&]() {
      resultEigen = solveWithEigen(N, b);
    });

    std::vector<double> resultGauss;
    double timeGauss = measureExecutionTime([&]() {
      resultGauss = solveWithGauss(N, b);
    });

    std::vector<double> resultQR;
    double timeQR = measureExecutionTime([&]() {
      resultQR = solveWithQR(N, b);
    });

    dataFile << N << " " << timeSherman << " " << timeEigen << " "
             << timeGauss << " " << timeQR << std::endl;

    std::cout << "Rozmiar macierzy: " << N << std::endl
              << "Czas Sherman-Morrison: " << timeSherman << " s\n"
              << "Czas Eigen: " << timeEigen << " s\n"
              << "Czas Gauss: " << timeGauss << " s\n"
              << "Czas QR: " << timeQR << " s\n";

    if (N == 120) {
      printVector(resultSherman, "Sherman-Morrison");
      printVector(resultEigen, "Eigen");
      printVector(resultGauss, "Gauss");
      printVector(resultQR, "QR");

      std::cout << "\nAnaliza błędów dla N = " << N << ":\n";
      std::cout << "Błąd między Sherman-Morrison a Eigen: "
                << calculateError(resultSherman, resultEigen) << std::endl;
      std::cout << "Błąd między Sherman-Morrison a Gauss: "
                << calculateError(resultSherman, resultGauss) << std::endl;
      std::cout << "Błąd między Sherman-Morrison a QR: "
                << calculateError(resultSherman, resultQR) << std::endl;
    }

    std::cout << "------------------------" << std::endl;
  }

  for (const auto& N : largerSizes) {
    std::vector<double> b(N, 2.0);
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 1.0));

    for (std::size_t i = 0; i < N; i++) {
      A[i][i] = 5.0;
      if (i < N - 1) {
        A[i][i + 1] = 3.0;
      }
    }

    double timeSherman = measureExecutionTime([&]() {
      solveShermanMorrison(A, b);
    });

    dataFileLarge << N << " " << timeSherman << std::endl;

    std::cout << "Rozmiar macierzy (duża): " << N << std::endl
              << "Czas Sherman-Morrison: " << timeSherman << " s\n"
              << "------------------------" << std::endl;
  }

  dataFile.close();
  dataFileLarge.close();

  createGnuplotScript("plot_script.gnu", "Porównanie czasów wykonania");
  createGnuplotScript("plot_script_large.gnu", "Czas wykonania Sherman-Morrison dla dużych macierzy", true);

  system("gnuplot plot_script.gnu");
  system("gnuplot plot_script_large.gnu");

  std::cout << "Wykresy zostały zapisane do plików 'timing_comparison.svg' i 'timing_comparison_large.svg'" << std::endl;

  return 0;
}
