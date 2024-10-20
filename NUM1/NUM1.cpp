#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

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

void process_data(double val, double initial_step, double min_step, double factor, const std::string& filename) {
  std::ofstream output(filename);
  output << "# Step_Size Forward_Error Central_Error\n";

  for (double step_size = initial_step; step_size >= min_step; step_size /= factor) {
    double exact_value = function_derivative(val);
    double forward_err = std::abs(forward_difference(val, step_size) - exact_value);
    double central_err = std::abs(central_difference(val, step_size) - exact_value);

    output << step_size << " " << forward_err << " " << central_err << "\n";
  }
  output.close();
}

void generate_plot(const std::string& filename) {
  FILE* pipe = popen("gnuplot -persist", "w");
  if (pipe) {
    fprintf(pipe, "set logscale xy\n");
    fprintf(pipe, "set title 'Error in Numerical Differentiation'\n");
    fprintf(pipe, "set xlabel 'Step Size (h)'\n");
    fprintf(pipe, "set ylabel 'Error'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "plot '%s' using 1:2 with linespoints title 'Forward Difference', '%s' using 1:3 with linespoints title 'Central Difference'\n", filename.c_str(), filename.c_str());
    fflush(pipe);
    pclose(pipe);
  } else {
    std::cerr << "Failed to initiate gnuplot." << std::endl;
  }
}

int main() {
  double x = 0.2;
  double h_start = 0.1;
  double h_end = 1e-10;
  double step = 1.05;

  std::string output_file = "NUM1_results.dat";
  process_data(x, h_start, h_end, step, output_file);
  generate_plot(output_file);

  return 0;
}
