#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

double f(double x) {
  return 1.0 / (1.0 + 10.0 * x * x);
}

double f1(double x) {
  return 1.0 / (1.0 + 10.0 * x * x);
}

double f3(double x) {
  return sin(3 * x) + x;
}

double f4(double x) {
  return sin(x) + 0.3 * sin(50 * x);
}

double f6(double x) {
  return (sin(2 * M_PI * x) + 0.2 * sin(100 * x)) / (1.3 * x + 2);
}

double f7(double x) {
  return exp(-x * x) * sin(10 * x);
}

std::vector<double> generateGrid(int n, double a, double b) {
  std::vector<double> X(n + 1);
  double step = (b - a) / n;
  for (int i = 0; i <= n; i++) {
    X[i] = a + i * step;
  }
  return X;
}

std::vector<double> evaluateFunction(const std::vector<double>& X, double (*f)(double)) {
  std::vector<double> Y(X.size());
  for (size_t i = 0; i < X.size(); i++) {
    Y[i] = f(X[i]);
  }
  return Y;
}

double lagrangeInterpolation(const std::vector<double>& X, const std::vector<double>& Y, double x) {
  double result = 0.0;
  int n = static_cast<int>(X.size());

  for (int i = 0; i < n; i++) {
    double Li = 1.0;
    for (int j = 0; j < n; j++) {
      if (j != i) {
        Li *= (x - X[j]) / (X[i] - X[j]);
      }
    }
    result += Y[i] * Li;
  }
  return result;
}

std::vector<double> computeSplineSecondDerivatives(const std::vector<double>& X, const std::vector<double>& Y) {
  int n = static_cast<int>(X.size()) - 1;
  std::vector<double> M(n + 1, 0.0);

  std::vector<double> alpha(n, 0.0);
  std::vector<double> beta(n, 0.0);
  std::vector<double> h(n);

  for (int i = 0; i < n; i++) {
    h[i] = X[i + 1] - X[i];
  }

  for (int i = 1; i < n; i++) {
    double hi = h[i];
    double him1 = h[i - 1];
    double yi = (Y[i + 1] - Y[i]) / hi - (Y[i] - Y[i - 1]) / him1;
    alpha[i] = 3.0 * yi;
  }

  std::vector<double> c(n + 1, 0.0);
  std::vector<double> d(n + 1, 0.0);

  c[0] = 0.0;
  d[0] = 0.0;

  for (int i = 1; i < n; i++) {
    double diag = 2.0 * (h[i - 1] + h[i]);
    double lower = h[i - 1];
    double upper = h[i];

    double tmp = (lower * c[i - 1] + diag);
    c[i] = -upper / tmp;
    d[i] = (alpha[i] - lower * d[i - 1]) / tmp;
  }

  for (int i = n - 1; i > 0; i--) {
    M[i] = d[i] + c[i] * M[i + 1];
  }

  return M;
}

double splineValue(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& M, double x) {
  int n = static_cast<int>(X.size()) - 1;

  int iLeft = 0;
  int iRight = n;
  while (iRight - iLeft > 1) {
    int mid = (iLeft + iRight) / 2;
    if (x < X[mid]) {
      iRight = mid;
    } else {
      iLeft = mid;
    }
  }
  int i = iLeft;

  double h = X[i + 1] - X[i];
  double A = (X[i + 1] - x) / h;
  double B = 1.0 - A;
  double C = (1.0 / 6.0) * (A * A * A - A) * h * h;
  double D = (1.0 / 6.0) * (B * B * B - B) * h * h;

  double S = A * Y[i] + B * Y[i + 1] + C * M[i] + D * M[i + 1];
  return S;
}

int main() {
  system("mkdir -p data");
  system("mkdir -p images");

  std::vector<int> nodes = {10, 20, 30};
  std::vector<std::pair<double (*)(double), std::string>> functions = {
      {f1, "runge"},
      {f3, "sin3x"},
      {f4, "sin_composite"},
      {f6, "rational_sin"},
      {f7, "gauss_sin"}};

  for (const auto& func : functions) {
    for (int n : nodes) {
      double (*f_ptr)(double) = func.first;
      std::vector<double> X = generateGrid(n, -1.0, 1.0);
      std::vector<double> Y = evaluateFunction(X, f_ptr);
      std::vector<double> M = computeSplineSecondDerivatives(X, Y);

      std::string filename = "data/interpolation_" + func.second + "_n" +
                             std::to_string(n) + ".dat";
      std::ofstream outFile(filename);
      outFile << std::fixed << std::setprecision(6);

      int plotPoints = 200;
      double step = 2.0 / (plotPoints - 1);
      double xx = -1.0;

      for (int k = 0; k < plotPoints; k++) {
        double fx = f_ptr(xx);
        double wx = lagrangeInterpolation(X, Y, xx);
        double sx = splineValue(X, Y, M, xx);

        outFile << xx << " " << fx << " " << wx << " "
                << std::fabs(fx - wx) << " " << sx << " "
                << std::fabs(fx - sx) << "\n";

        xx += step;
      }
      outFile.close();
    }
  }

  std::ofstream plotScript("plot_functions.gnu");
  plotScript << "set terminal svg enhanced size 1200,800\n";
  plotScript << "set grid\n";
  plotScript << "set key outside right spacing 2\n\n";

  for (const auto& func : functions) {
    for (int n : nodes) {
      plotScript << "set output 'images/interpolation_" << func.second << "_n" << n << ".svg'\n";

      if (func.second == "sin3x" && n == 20) {
        plotScript << "set xrange [-2:2]\n";
        plotScript << "set yrange [-2:2]\n";
      } else {
        plotScript << "set xrange [-1:1]\n";
        plotScript << "set autoscale y\n";
      }

      plotScript << "set xlabel 'x'\n";
      plotScript << "set ylabel 'y'\n";

      std::string title;
      if (func.second == "runge") {
        title = "f(x) = 1/(1+10x^2)";
      } else if (func.second == "sin3x") {
        title = "f(x) = sin(3x) + x";
      } else if (func.second == "sin_composite") {
        title = "f(x) = sin(x) + 0.3sin(50x)";
      } else if (func.second == "gauss_sin") {
        title = "f(x) = e^{-x^2}sin(10x)";
      } else {
        title = "f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2)";
      }

      plotScript << "set title '" << title << ", n = " << n << " węzłów'\n";

      plotScript << "plot 'data/interpolation_" << func.second << "_n" << n
                 << ".dat' using 1:2 title 'Funkcja oryginalna' w lines lc rgb '#0060ad', \\\n"
                 << "     '' using 1:3 title 'Interpolacja Lagrange''a' w lines lc rgb '#dd181f', \\\n"
                 << "     '' using 1:5 title 'Spline kubiczny' w lines lc rgb '#00A000'\n\n";
    }
  }
  plotScript.close();

  std::ofstream errorScript("plot_errors.gnu");
  errorScript << "set terminal svg enhanced size 1200,800\n";
  errorScript << "set grid\n";
  errorScript << "set key outside right spacing 2\n";
  errorScript << "set logscale y\n";
  errorScript << "set format y '10^{%L}'\n";
  errorScript << "set xrange [-1:1]\n";
  errorScript << "set yrange [1e-5:1e0]\n";
  errorScript << "set lmargin 12\n";
  errorScript << "set rmargin 6\n";
  errorScript << "set ylabel offset -2,0\n";
  errorScript << "set ytics offset 0,0\n\n";

  for (const auto& func : functions) {
    for (int n : nodes) {
      errorScript << "set output 'images/error_" << func.second << "_n" << n << ".svg'\n";

      std::string title;
      if (func.second == "runge") {
        title = "f(x) = 1/(1+10x^2)";
      } else if (func.second == "sin3x") {
        title = "f(x) = sin(3x) + x";
      } else if (func.second == "sin_composite") {
        title = "f(x) = sin(x) + 0.3sin(50x)";
      } else if (func.second == "gauss_sin") {
        title = "f(x) = e^{-x^2}sin(10x)";
      } else {
        title = "f(x) = (sin(2πx) + 0.2sin(100x))/(1.3x + 2)";
      }

      errorScript << "set title 'Błędy interpolacji funkcji " << title << "\\n = " << n << "'\n";
      errorScript << "set xlabel 'x'\n";
      errorScript << "set ylabel 'Błąd' offset -3\n";

      errorScript << "plot 'data/interpolation_" << func.second << "_n" << n
                  << ".dat' using 1:4 title 'Błąd Lagrange''a' w lines lc rgb '#dd181f', \\\n"
                  << "     '' using 1:6 title 'Błąd splinu' w lines lc rgb '#00A000'\n\n";
    }
  }
  errorScript << "unset logscale y\n";
  errorScript.close();

  system("gnuplot plot_functions.gnu");
  system("gnuplot plot_errors.gnu");

  return 0;
}
