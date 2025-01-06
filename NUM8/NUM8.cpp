#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "helper.h"

double safeLogSquared(double x) {
  double val = x * x;
  if (val < 1e-15) {
    val = 1e-15;
  }
  return std::log(val);
}

double createFunction0(double x) {
  return std::sin(x * std::exp(x));
}

double createFunction1(double x) {
  double logVal = safeLogSquared(x);
  return std::cos(logVal * std::tanh(x));
}

double createFunction2(double x) {
  return std::exp(x * x * std::sin(x * x));
}

double createFunction3(double x) {
  return std::exp(std::sin(5 * x) * std::sinh(x));
}

double createFunction4(double x) {
  return std::sin(std::cos(std::exp(x)));
}

double createFunction5(double x) {
  return std::cos(1 / (1 + std::abs(x)) * std::exp(x * x));
}

double createFunctionMap(double x, int i) {
  switch (i) {
    case 0:
      return createFunction0(x);
    case 1:
      return createFunction1(x);
    case 2:
      return createFunction2(x);
    case 3:
      return createFunction3(x);
    case 4:
      return createFunction4(x);
    case 5:
      return createFunction5(x);
    default:
      return 0.0;
  }
}

double function(double x, const std::vector<double>& a) {
  double result = 0.0;
  int M = static_cast<int>(a.size());
  for (int i = 0; i < M; i++) {
    result += a[i] * createFunctionMap(x, i);
  }
  return result;
}

void functionApproximation(const int N, const std::vector<double>& a, double sigma) {
  std::cout << "Zaburzenie (sigma): " << sigma << std::endl;
  std::cout << "Liczba punktów (N): " << N << std::endl;

  int M = static_cast<int>(a.size());

  std::random_device rd;
  std::mt19937 generator(rd());
  std::normal_distribution<double> normalDist(0.0, sigma);

  std::vector<double> xPoints;
  xPoints.reserve(N);
  std::vector<double> yPoints;
  yPoints.reserve(N);

  for (int i = 0; i < N; i++) {
    double t = i / static_cast<double>(N - 1);
    double x = -1.0 + 2.0 * t;
    xPoints.push_back(x);

    double noise = normalDist(generator);

    yPoints.push_back(function(x, a) + noise);
  }

  std::ofstream dataFile("points.dat");
  for (size_t i = 0; i < xPoints.size(); i++) {
    dataFile << xPoints[i] << " " << yPoints[i] << "\n";
  }
  dataFile.close();

  std::vector<std::vector<double>> matrix(N, std::vector<double>(M));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      matrix[i][j] = createFunctionMap(xPoints[i], j);
    }
  }

  Eigen::MatrixXd eigenMatrix = convertMatrixToEigen(matrix);
  Eigen::VectorXd eigenY = convertVectorToEigen(yPoints);

  Eigen::VectorXd eigenCoefficients =
      eigenMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(eigenY);

  std::cout << "Wyznaczone wspolczynniki aproksymacji:" << std::endl;
  printEigenVector(eigenCoefficients);

  std::vector<double> approxCoeff = convertEigenToVector(eigenCoefficients);

  std::ofstream approxFile("approximation.dat");
  const int approxPoints = 1000;
  for (int i = 0; i < approxPoints; i++) {
    double t = i / static_cast<double>(approxPoints - 1);
    double x = -1.0 + 2.0 * t;
    approxFile << x << " " << function(x, approxCoeff) << "\n";
  }
  approxFile.close();

  std::ofstream originalFile("original.dat");
  const int originalPoints = 1000;
  for (int i = 0; i < originalPoints; i++) {
    double t = i / static_cast<double>(originalPoints - 1);
    double x = -1.0 + 2.0 * t;
    originalFile << x << " " << function(x, a) << "\n";
  }
  originalFile.close();

  std::string plotFileName = "plot_N_" + std::to_string(N) + "_Sigma_" +
                             std::to_string(sigma) + ".gnu";
  std::ofstream plotScript(plotFileName);

  plotScript << "set terminal svg enhanced size 800,600\n";
  plotScript << "set output 'plot_N_" << N << "_Sigma_" << sigma << ".svg'\n";
  plotScript << "set title 'Aproksymacja (N=" << N << ", σ=" << sigma << ")'\n";
  plotScript << "set xlabel 'x'\n";
  plotScript << "set ylabel 'y'\n";
  plotScript << "set grid\n";
  plotScript << "set xrange [-1:1]\n";
  plotScript << "plot 'points.dat' using 1:2 with points pt 7 ps 0.5 title 'Punkty pomiarowe', \\\n";
  plotScript << "     'original.dat' using 1:2 with lines lw 2 title 'Funkcja F(x)', \\\n";
  plotScript << "     'approximation.dat' using 1:2 with lines lw 2 dt 2 title 'Aproksymacja'\n";
  plotScript.close();

  system(("gnuplot " + plotFileName).c_str());
}

int main() {
  std::vector<double> coefficients = {-0.5, 1.0, -0.5, 1.5, -0.5, -2.5};

  functionApproximation(10, coefficients, 0.1);
  functionApproximation(25, coefficients, 0.1);
  functionApproximation(70, coefficients, 0.1);
  functionApproximation(100, coefficients, 0.1);

  functionApproximation(10, coefficients, 0.5);
  functionApproximation(25, coefficients, 0.5);
  functionApproximation(70, coefficients, 0.5);
  functionApproximation(100, coefficients, 0.5);

  functionApproximation(10, coefficients, 1.5);
  functionApproximation(25, coefficients, 1.5);
  functionApproximation(70, coefficients, 1.5);
  functionApproximation(100, coefficients, 1.5);
  return 0;
}
