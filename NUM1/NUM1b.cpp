#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

float function(float x) {
  return sin(pow(x, 3));  // f(x) = sin(x^3)
}

double function(double x) {
  return sin(pow(x, 3));  // f(x) = sin(x^3)
}

float function_derivative(float x) {
  return 3 * pow(x, 2) * cos(pow(x, 3));  // f'(x) = 3x^2cos(x^3)
}

double function_derivative(double x) {
  return 3 * pow(x, 2) * cos(pow(x, 3));  // f'(x) = 3x^2cos(x^3)
}

float central_difference(float x, float h) {
  return (function(x + h) - function(x - h)) / (2 * h);
}

double central_difference(double x, double h) {
  return (function(x + h) - function(x - h)) / (2 * h);
}

void plot_data(const std::vector<float>& step_sizes_float, const std::vector<float>& central_arr_float,
               const std::vector<double>& step_sizes_double, const std::vector<double>& central_arr_double) {
  FILE* pipe = popen("gnuplot -persist", "w");
  if (pipe) {
    fprintf(pipe, "set logscale xy\n");
    fprintf(pipe, "set title 'Error in Central Difference (float vs double)'\n");
    fprintf(pipe, "set xlabel 'Step Size (h)'\n");
    fprintf(pipe, "set ylabel 'Error'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "plot '-' using 1:2 with lines title 'Central Difference (float)', '-' using 1:2 with lines title 'Central Difference (double)'\n");

    for (size_t i = 0; i < step_sizes_float.size(); ++i) {
      fprintf(pipe, "%g %g\n", step_sizes_float[i], central_arr_float[i]);
    }
    fprintf(pipe, "e\n");

    for (size_t i = 0; i < step_sizes_double.size(); ++i) {
      fprintf(pipe, "%g %g\n", step_sizes_double[i], central_arr_double[i]);
    }
    fprintf(pipe, "e\n");

    fflush(pipe);
    pclose(pipe);
  } else {
    std::cerr << "Failed to initiate gnuplot." << std::endl;
  }
}

int main() {
  float x_float = 0.2;
  double x_double = 0.2;
  float h_start = 0.1;
  float h_end = 1e-18;
  float factor = 1.1;

  std::vector<float> central_arr_float;
  std::vector<float> step_sizes_float;

  for (float h = h_start; h >= h_end; h /= factor) {
    float exact_value_float = function_derivative(x_float);
    float central_err_float = std::abs(central_difference(x_float, h) - exact_value_float);

    central_arr_float.push_back(central_err_float);
    step_sizes_float.push_back(h);
  }

  std::vector<double> central_arr_double;
  std::vector<double> step_sizes_double;

  for (double h = h_start; h >= h_end; h /= factor) {
    double exact_value_double = function_derivative(x_double);
    double central_err_double = std::abs(central_difference(x_double, h) - exact_value_double);

    central_arr_double.push_back(central_err_double);
    step_sizes_double.push_back(h);
  }

  plot_data(step_sizes_float, central_arr_float, step_sizes_double, central_arr_double);

  return 0;
}
