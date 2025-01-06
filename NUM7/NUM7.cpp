#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "helper.h"

double y(double x) {
  return 1.0 / (1.0 + 10.0 * x * x);
}

double f1(double x) {
  return 12 * std::pow(x, 5) - 5 * std::pow(x, 3);
}

double f2(double x) {
  return std::sin(2 * M_E * std::pow(x, 2)) / (2 + M_PI * x * x);
}

double f3(double x) {
  return std::exp(-x * x);
}

double lagrangeInterpolation(double x, std::vector<double> xPoints, std::vector<double> yPoints, int n) {
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    double term = yPoints[i];
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
        term = term * (x - xPoints[j]) / (xPoints[i] - xPoints[j]);
      }
    }
    result += term;
  }
  return result;
}

std::vector<std::vector<double>> cubicSpline(std::vector<double> &xPoints, std::vector<double> &yPoints, int n) {
  double h = xPoints[1] - xPoints[0];

  std::vector<double> epsilons(n);
  for (size_t i = 1; i < n - 1; i++) {
    epsilons[i] = (6 * (yPoints[i - 1] - 2 * yPoints[i] + yPoints[i + 1])) / (h * h);
  }
  epsilons[0] = 0;
  epsilons[n - 1] = 0;

  std::vector<double> diagonal(n, 4.0);
  std::vector<double> subDiagonals(n - 1, 1.0);
  std::vector<double> superDiagonals(n - 1, 1.0);
  diagonal[0] = 1.0;
  diagonal[n - 1] = 1.0;

  for (int i = 1; i < n; i++) {
    double factor = subDiagonals[i - 1] / diagonal[i - 1];
    diagonal[i] -= factor * superDiagonals[i - 1];
    epsilons[i] -= factor * epsilons[i - 1];
  }

  std::vector<double> solution(n);
  solution[n - 1] = epsilons[n - 1] / diagonal[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    solution[i] = (epsilons[i] - superDiagonals[i] * solution[i + 1]) / diagonal[i];
  }

  std::vector<std::vector<double>> coefficients(4, std::vector<double>(n - 1, 0.0));
  for (int i = 0; i < n - 1; i++) {
    coefficients[0][i] = (solution[i + 1] - solution[i]) / (6.0 * h);
    coefficients[1][i] = solution[i] / 2.0;
    coefficients[2][i] = (yPoints[i + 1] - yPoints[i]) / h - (h / 6.0) * (2.0 * solution[i] + solution[i + 1]);
    coefficients[3][i] = yPoints[i];
  }

  return coefficients;
}

void task(std::string name, int N) {
  std::map<std::string, std::function<double(double)>> functions;
  functions["f1"] = f1;
  functions["f2"] = f2;
  functions["f3"] = f3;

  auto func = functions[name];

  std::vector<double> xPoints;
  xPoints.reserve(N);
  std::vector<double> yPoints;
  yPoints.reserve(N);

  for (int i = 0; i < N; i++) {
    double t = i / static_cast<double>(N - 1);
    double x = -1.0 + 2.0 * t;
    xPoints.push_back(x);
    yPoints.push_back(func(x));
  }

  for (size_t i = 0; i < N; i++) {
    double x = xPoints[i];
    double value = lagrangeInterpolation(x, xPoints, yPoints, N);
    std::cout << name << "(" << x << ") = " << value << std::endl;
  }
}

int main() {
  task("f1", 10);
  task("f2", 10);
  task("f2", 10);
  return 0;
}