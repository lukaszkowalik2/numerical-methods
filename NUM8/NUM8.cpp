#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "helper.h"

double createFunction0(double x) {
  return std::sin(x * std::exp(x));
}

double createFunction1(double x) {
  return std::cos(std::log(x * x) * std::tanh(x));
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

double distribution(double value, double sigma) {
  return (1 / (std::sqrt(2 * M_PI * sigma * sigma))) * std::exp(-(value * value) / (2 * sigma * sigma));
}

double function(double x, std::vector<double> a) {
  double result = 0.0;
  int M = a.size();
  for (int i = 0; i < M; i++) {
    result += a[i] * createFunctionMap(x, i);
  }
  return result;
}

void functionApproximation(const int N, std::vector<double> a, double sigma) {
  int M = a.size();

  std::vector<double> xPoints;
  std::vector<double> yPoints;
  for (int i = 0; i < N; i++) {
    double t = i / static_cast<double>(N - 1);
    xPoints.push_back(-1.0 + 2.0 * t);
    double noise = distribution(xPoints[i], sigma);
    yPoints.push_back(function(xPoints[i], a) + noise);
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
  Eigen::VectorXd eigenCoefficients = eigenMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(convertVectorToEigen(yPoints));
  printEigenVector(eigenCoefficients);

  std::ofstream approxFile("approximation.dat");
  std::vector<double> approxCoeff = convertEigenToVector(eigenCoefficients);
  for (int i = 0; i < N; i++) {
    approxFile << xPoints[i] << " " << function(xPoints[i], approxCoeff) << "\n";
  }
  approxFile.close();

  std::ofstream originalFile("original.dat");
  for (int i = 0; i < N; i++) {
    originalFile << xPoints[i] << " " << function(xPoints[i], a) << "\n";
  }
  originalFile.close();

  system("gnuplot plot.gnu");
}

int main() {
  std::vector<double> coefficients = {-0.5, 1.0, -0.5, 1.5, -0.5, -2.5};
  functionApproximation(100, coefficients, 0.1);
  return 0;
}
