#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

double y(double x) {
  return 1.0 / (1.0 + 10.0 * x * x);
}

double f1(double x) {
  return 12.0 * std::pow(x, 5) - 5.0 * std::pow(x, 3);
}

double f2(double x) {
  return std::sin(2 * M_E * x * x) / (2 + M_PI * x * x);
}

double f3(double x) {
  return std::exp(-x * x);
}

double lagrangeInterpolation(double x, std::vector<double> &xPoints, std::vector<double> &yPoints) {
  const size_t n = xPoints.size();
  double result = 0.0;
  for (size_t i = 0; i < n; i++) {
    double term = yPoints[i];
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
        term *= (x - xPoints[j]) / (xPoints[i] - xPoints[j]);
      }
    }
    result += term;
  }
  return result;
}

std::vector<std::vector<double>> cubicSpline(std::vector<double> &xPoints, std::vector<double> &yPoints, int n) {
  double h = xPoints[1] - xPoints[0];

  std::vector<double> epsilons(n, 0.0);
  for (int i = 1; i < n - 1; i++) {
    epsilons[i] = 6.0 * (yPoints[i - 1] - 2.0 * yPoints[i] + yPoints[i + 1]) / (h * h);
  }

  std::vector<double> diagonal(n, 4.0);
  std::vector<double> subDiagonal(n - 1, 1.0);
  std::vector<double> superDiagonal(n - 1, 1.0);
  diagonal[0] = 1.0;
  diagonal[n - 1] = 1.0;

  for (int i = 1; i < n; i++) {
    double factor = subDiagonal[i - 1] / diagonal[i - 1];
    diagonal[i] -= factor * superDiagonal[i - 1];
    epsilons[i] -= factor * epsilons[i - 1];
  }

  std::vector<double> M(n, 0.0);
  M[n - 1] = epsilons[n - 1] / diagonal[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    M[i] = (epsilons[i] - superDiagonal[i] * M[i + 1]) / diagonal[i];
  }

  std::vector<std::vector<double>> coefficients(4, std::vector<double>(n - 1, 0.0));
  for (int i = 0; i < n - 1; i++) {
    coefficients[0][i] = (M[i + 1] - M[i]) / (6.0 * h);
    coefficients[1][i] = M[i] / 2.0;
    coefficients[2][i] = (yPoints[i + 1] - yPoints[i]) / h - (h / 6.0) * (2.0 * M[i] + M[i + 1]);
    coefficients[3][i] = yPoints[i];
  }
  return coefficients;
}

double calculateSpline(double x, std::vector<double> &xPoints, std::vector<std::vector<double>> &coefficients, int n) {
  int iInterval = n - 2;
  for (int i = 0; i < n - 1; i++) {
    if (x >= xPoints[i] && x <= xPoints[i + 1]) {
      iInterval = i;
      break;
    }
  }

  double a = coefficients[0][iInterval];
  double b = coefficients[1][iInterval];
  double c = coefficients[2][iInterval];
  double d = coefficients[3][iInterval];
  double dx = x - xPoints[iInterval];
  return a * std::pow(dx, 3) + b * std::pow(dx, 2) + c * dx + d;
}

void task(const std::string &funcName, int N) {
  // Create charts directory if it doesn't exist
  std::filesystem::create_directory("charts");

  static std::map<std::string, std::function<double(double)>> functions = {
      {"y", y},
      {"f1", f1},
      {"f2", f2},
      {"f3", f3}};

  if (functions.find(funcName) == functions.end()) {
    std::cerr << "Unknown function name: " << funcName << std::endl;
    return;
  }

  auto func = functions[funcName];
  std::vector<double> xPoints(N), yPoints(N);
  for (int i = 0; i < N; i++) {
    double t = i / double(N - 1);
    double x = -1.0 + 2.0 * t;
    xPoints[i] = x;
    yPoints[i] = func(x);
  }

  auto splineCoefficients = cubicSpline(xPoints, yPoints, N);

  int M = 200;
  std::string csvFile = "charts/results_" + funcName + "_" + std::to_string(N) + ".csv";
  std::ofstream fout(csvFile);
  fout << "# x, f(x), Lagrange, errLag, Spline, errSpl\n";

  std::vector<double> X(M), F(M), FLag(M), FSpl(M);

  for (int i = 0; i < M; i++) {
    double t = i / double(M - 1);
    double x = -1.0 + 2.0 * t;
    double fExact = func(x);
    double fLag = lagrangeInterpolation(x, xPoints, yPoints);
    double fSpl = calculateSpline(x, xPoints, splineCoefficients, N);
    double errLag = std::fabs(fExact - fLag);
    double errSpl = std::fabs(fExact - fSpl);
    fout << x << ","
         << fExact << ","
         << fLag << ","
         << errLag << ","
         << fSpl << ","
         << errSpl << "\n";
    X[i] = x;
    F[i] = fExact;
    FLag[i] = fLag;
    FSpl[i] = fSpl;
  }
  fout.close();

  std::cout << "Data saved to " << csvFile << std::endl;

  std::string gpFile = "charts/plot_" + funcName + "_" + std::to_string(N) + ".gp";
  std::string svgFile = "charts/plot_" + funcName + "_" + std::to_string(N) + ".svg";
  {
    std::ofstream gp(gpFile);
    gp << "set terminal svg size 1200,800 enhanced font 'Arial,12'\n";
    gp << "set output \"" << svgFile << "\"\n";
    gp << "set title \"Interpolation of " << funcName
       << " with N=" << N << "\" font 'Arial,14'\n";
    gp << "set xlabel \"x\" font 'Arial,12'\n";
    gp << "set ylabel \"f(x)\" font 'Arial,12'\n";
    gp << "set grid lw 1\n";
    gp << "set key outside right\n";
    gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2\n";
    gp << "set style line 2 lc rgb '#dd181f' lt 1 lw 2\n";
    gp << "set style line 3 lc rgb '#00cc00' lt 1 lw 2\n";
    gp << "plot \"" << csvFile << "\" using 1:2 with lines ls 1 title 'Exact', \\\n"
       << "     \"" << csvFile << "\" using 1:3 with lines ls 2 title 'Lagrange', \\\n"
       << "     \"" << csvFile << "\" using 1:5 with lines ls 3 title 'Spline'\n";
    gp << "unset output\n";
    gp << "exit\n";
  }

  std::string cmd = "gnuplot " + gpFile;
  int ret = system(cmd.c_str());
  if (ret != 0) {
    std::cerr << "Warning: gnuplot command failed or not found.\n";
  } else {
    std::cout << "SVG plot generated: " << svgFile << std::endl;
    std::filesystem::remove(gpFile);
  }
}

int main() {
  task("f1", 10);
  task("f2", 10);
  task("f3", 10);

  return 0;
}
