#include <stdio.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double function(double x) {
  return sin(pow(x, 3));  // f(x) = sin(x^3)
}

double function_derivative(double x) {
  return 3 * pow(x, 2) * cos(pow(x, 3));  // f'(x) = 3x^2cos(x^3)
}

double forward_difference(double x, double h) {
  return (function(x + h) - function(x)) / h;
}

double central_difference(double x, double h) {
  return (function(x + h) - function(x - h)) / (2 * h);
}

void process_data_double(double val, double initial_step, double min_step, double factor, const std::string& filename) {
  std::vector<double> forward_arr;
  std::vector<double> central_arr;
  std::ofstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  for (double step_size = initial_step; step_size >= min_step; step_size /= factor) {
    double exact_value = function_derivative(val);
    double forward_err = std::abs(forward_difference(val, step_size) - exact_value);
    double central_err = std::abs(central_difference(val, step_size) - exact_value);

    forward_arr.push_back(forward_err);
    central_arr.push_back(central_err);

    file << step_size << " " << forward_err << " " << central_err << std::endl;
  }

  file.close();
}

void process_data_float(float val, float initial_step, float min_step, float factor, const std::string& filename) {
  std::vector<float> forward_arr;
  std::vector<float> central_arr;
  std::ofstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  for (float step_size = initial_step; step_size >= min_step; step_size /= factor) {
    float exact_value = function_derivative(val);
    float forward_err = std::abs(forward_difference(val, step_size) - exact_value);
    float central_err = std::abs(central_difference(val, step_size) - exact_value);

    forward_arr.push_back(forward_err);
    central_arr.push_back(central_err);

    file << step_size << " " << forward_err << " " << central_err << std::endl;
  }

  file.close();
}

void generate_plot(const std::string& filename_double, const std::string& filename_float) {
  FILE* pipe = popen("gnuplot -persist", "w");
  if (pipe) {
    fprintf(pipe, "set logscale xy\n");
    fprintf(pipe, "set multiplot layout 2,1 title 'Error in Numerical Differentiation'\n");

    fprintf(pipe, "set title 'Double Precision'\n");
    fprintf(pipe, "set xlabel 'Step Size (h)'\n");
    fprintf(pipe, "set ylabel 'Error'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "plot '%s' using 1:2 with linespoints title 'Forward Difference (double)', '%s' using 1:3 with linespoints title 'Central Difference (double)'\n", filename_double.c_str(), filename_double.c_str());

    fprintf(pipe, "set title 'Float Precision'\n");
    fprintf(pipe, "set xlabel 'Step Size (h)'\n");
    fprintf(pipe, "set ylabel 'Error'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "plot '%s' using 1:2 with linespoints title 'Forward Difference (float)', '%s' using 1:3 with linespoints title 'Central Difference (float)'\n", filename_float.c_str(), filename_float.c_str());

    fprintf(pipe, "unset multiplot\n");
    fflush(pipe);
    pclose(pipe);
  } else {
    std::cerr << "Failed to initiate gnuplot." << std::endl;
  }
}

int main() {
  double x_double = 0.2;
  double h_start_double = 0.1;
  double h_end_double = 1e-10;
  double step_double = 1.1;
  std::string filename_double = "data_double.txt";

  float x_float = (float)x_double;
  float h_start_float = (float)h_start_double;
  float h_end_float = (float)h_end_double;
  float step_float = (float)step_double;
  std::string filename_float = "data_float.txt";

  process_data_double(x_double, h_start_double, h_end_double, step_double, filename_double);
  process_data_float(x_float, h_start_float, h_end_float, step_float, filename_float);
  generate_plot(filename_double, filename_float);

  return 0;
}