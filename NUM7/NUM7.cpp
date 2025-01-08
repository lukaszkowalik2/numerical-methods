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
  return -std::log(std::abs(std::sin(M_PI * x) + 0.1)) / 4.0;
}

double lagrangeInterpolation(double x, std::vector<double>& xPoints, std::vector<double>& yPoints) {
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

std::vector<double> solveTridiagonalSystem(const std::vector<double>& diagonal,
                                           const std::vector<double>& subDiagonal,
                                           const std::vector<double>& superDiagonal,
                                           const std::vector<double>& rhs) {
  int n = diagonal.size();
  std::vector<double> solution(n);
  std::vector<double> temp_diag = diagonal;
  std::vector<double> temp_rhs = rhs;

  for (int i = 1; i < n; i++) {
    double factor = subDiagonal[i - 1] / temp_diag[i - 1];
    temp_diag[i] -= factor * superDiagonal[i - 1];
    temp_rhs[i] -= factor * temp_rhs[i - 1];
  }

  solution[n - 1] = temp_rhs[n - 1] / temp_diag[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    solution[i] = (temp_rhs[i] - superDiagonal[i] * solution[i + 1]) / temp_diag[i];
  }
  return solution;
}

std::map<char, std::vector<double>> cubicSpline(std::vector<double>& xPoints,
                                                std::vector<double>& yPoints,
                                                int n) {
  std::map<char, std::vector<double>> coefficients;
  double h = xPoints[1] - xPoints[0];

  std::vector<double> diagonal(n, 4.0);
  std::vector<double> subDiagonal(n - 1, 1.0);
  std::vector<double> superDiagonal(n - 1, 1.0);
  std::vector<double> rhs(n, 0.0);

  diagonal[0] = diagonal[n - 1] = 1.0;
  for (int i = 1; i < n - 1; i++) {
    rhs[i] = 6.0 * (yPoints[i - 1] - 2.0 * yPoints[i] + yPoints[i + 1]) / (h * h);
  }

  std::vector<double> M = solveTridiagonalSystem(diagonal, subDiagonal, superDiagonal, rhs);

  coefficients['a'] = std::vector<double>(n - 1);
  coefficients['b'] = std::vector<double>(n - 1);
  coefficients['c'] = std::vector<double>(n - 1);
  coefficients['d'] = std::vector<double>(n - 1);

  for (int i = 0; i < n - 1; i++) {
    coefficients['a'][i] = (M[i + 1] - M[i]) / (6.0 * h);
    coefficients['b'][i] = M[i] / 2.0;
    coefficients['c'][i] = (yPoints[i + 1] - yPoints[i]) / h - (h / 6.0) * (2.0 * M[i] + M[i + 1]);
    coefficients['d'][i] = yPoints[i];
  }

  return coefficients;
}

double calculateSpline(double x, std::vector<double>& xPoints,
                       std::map<char, std::vector<double>>& coef, int n) {
  int iInterval = n - 2;
  for (int i = 0; i < n - 1; i++) {
    if (x >= xPoints[i] && x <= xPoints[i + 1]) {
      iInterval = i;
      break;
    }
  }

  double dx = x - xPoints[iInterval];
  return coef['a'][iInterval] * std::pow(dx, 3) +
         coef['b'][iInterval] * std::pow(dx, 2) +
         coef['c'][iInterval] * dx +
         coef['d'][iInterval];
}

void task(const std::string& funcName, int N) {
  std::string outputDir = "output";
  std::filesystem::create_directories(outputDir);

  static std::map<std::string, std::function<double(double)>> functions = {
      {"y", y}, {"f1", f1}, {"f2", f2}, {"f3", f3}};

  if (functions.find(funcName) == functions.end()) {
    std::cerr << "Unknown function name: " << funcName << std::endl;
    return;
  }

  auto func = functions[funcName];
  std::vector<double> xPoints(N), yPoints(N);
  for (int i = 0; i < N; i++) {
    double t = i / double(N - 1);
    xPoints[i] = -1.0 + 2.0 * t;
    yPoints[i] = func(xPoints[i]);
  }

  auto splineCoeff = cubicSpline(xPoints, yPoints, N);

  std::ofstream lagrangeFile(outputDir + "/lagrange_" + funcName + "_" + std::to_string(N) + ".dat");
  std::ofstream splineFile(outputDir + "/spline_" + funcName + "_" + std::to_string(N) + ".dat");
  std::ofstream exactFile(outputDir + "/exact_" + funcName + "_" + std::to_string(N) + ".dat");
  std::ofstream errorFile(outputDir + "/errors_" + funcName + "_" + std::to_string(N) + ".dat");

  for (double x = -1.0; x <= 1.0; x += 0.001) {
    double exact = func(x);
    double lagrange = lagrangeInterpolation(x, xPoints, yPoints);
    double spline = calculateSpline(x, xPoints, splineCoeff, N);

    exactFile << std::scientific << x << " " << exact << "\n";
    lagrangeFile << std::scientific << x << " " << lagrange << "\n";
    splineFile << std::scientific << x << " " << spline << "\n";
    errorFile << std::scientific << x << " "
              << std::abs(lagrange - exact) << " "
              << std::abs(spline - exact) << "\n";
  }

  std::ofstream pointsFile(outputDir + "/points_" + funcName + "_" + std::to_string(N) + ".dat");
  for (int i = 0; i < N; i++) {
    pointsFile << std::scientific << xPoints[i] << " " << yPoints[i] << "\n";
  }

  lagrangeFile.close();
  splineFile.close();
  exactFile.close();
  errorFile.close();
  pointsFile.close();
  std::ofstream plotScript(outputDir + "/plot_" + funcName + "_" + std::to_string(N) + ".gnu");
  plotScript << "set terminal svg enhanced size 800,600\n";
  plotScript << "set output 'interpolation_" << funcName << "_" << N << ".svg'\n";
  plotScript << "set title 'Interpolation for " << funcName << " (N=" << N << ")'\n";
  plotScript << "set xlabel 'x'\n";
  plotScript << "set ylabel 'y'\n";
  plotScript << "set grid\n";

  // if ((funcName == "y" && N == 30) || (funcName == "f3" && (N == 20 || N == 30))) {
  //   plotScript << "set yrange [-2:2]\n";
  // }

  plotScript << "plot 'exact_" << funcName << "_" << N << ".dat' w l title 'Exact', \\\n"
             << "     'lagrange_" << funcName << "_" << N << ".dat' w l title 'Lagrange', \\\n"
             << "     'spline_" << funcName << "_" << N << ".dat' w l title 'Spline', \\\n"
             << "     'points_" << funcName << "_" << N << ".dat' w p pt 7 title 'Points'\n";
  plotScript.close();

  std::ofstream errorPlotScript(outputDir + "/plot_error_" + funcName + "_" + std::to_string(N) + ".gnu");
  errorPlotScript << "set terminal svg enhanced size 800,600\n";
  errorPlotScript << "set output 'errors_" << funcName << "_" << N << ".svg'\n";
  errorPlotScript << "set title 'Interpolation Errors for " << funcName << " (N=" << N << ")'\n";
  errorPlotScript << "set xlabel 'x'\n";
  errorPlotScript << "set ylabel 'Error'\n";
  errorPlotScript << "set grid\n";
  errorPlotScript << "set logscale y\n";
  errorPlotScript << "plot 'errors_" << funcName << "_" << N << ".dat' using 1:2 w l title 'Lagrange Error', \\\n"
                  << "     'errors_" << funcName << "_" << N << ".dat' using 1:3 w l title 'Spline Error'\n";
  errorPlotScript.close();

  std::string currentDir = std::filesystem::current_path().string();
  system(("cd " + currentDir + "/" + outputDir + " && gnuplot plot_" + funcName + "_" + std::to_string(N) + ".gnu").c_str());
  system(("cd " + currentDir + "/" + outputDir + " && gnuplot plot_error_" + funcName + "_" + std::to_string(N) + ".gnu").c_str());
}

int main() {
  task("y", 10);
  task("y", 20);
  task("y", 30);

  task("f1", 10);
  task("f1", 20);
  task("f1", 30);

  task("f2", 10);
  task("f2", 20);
  task("f2", 30);

  task("f3", 10);
  task("f3", 20);
  task("f3", 30);
  return 0;
}
