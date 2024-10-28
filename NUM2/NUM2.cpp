#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <random>

void solveAndAnalyze(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double perturbationNorm) {
  Eigen::VectorXd y = A.colPivHouseholderQr().solve(b);
  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Solution for A*y = b without perturbation:\n"
            << y.transpose() << "\n\n";

  Eigen::VectorXd delta_b = Eigen::VectorXd::Random(b.size()).normalized() * perturbationNorm;
  Eigen::VectorXd b_perturbed = b + delta_b;

  Eigen::VectorXd y_perturbed = A.colPivHouseholderQr().solve(b_perturbed);
  std::cout << "Solution for A*y = b + Δb (perturbed vector):\n"
            << y_perturbed.transpose() << "\n\n";

  Eigen::VectorXd difference = y_perturbed - y;
  std::cout << "Difference between perturbed and original solution:\n"
            << difference.transpose() << "\n";
  std::cout << "Norm of the difference: " << std::scientific << difference.norm() << "\n\n";
}

int main() {
  Eigen::MatrixXd A1(5, 5);
  A1 << 5.8267103432, 1.0419816676, 0.4517861296, -0.2246976350, 0.7150286064,
      1.0419816676, 5.8150823499, -0.8642832971, 0.6610711416, -0.3874139415,
      0.4517861296, -0.8642832971, 1.5136472691, -0.8512078774, 0.6771688230,
      -0.2246976350, 0.6610711416, -0.8512078774, 5.3014166511, 0.5228116055,
      0.7150286064, -0.3874139415, 0.6771688230, 0.5228116055, 3.5431433879;

  Eigen::MatrixXd A2(5, 5);
  A2 << 5.4763986379, 1.6846933459, 0.3136661779, -1.0597154562, 0.0083249547,
      1.6846933459, 4.6359087874, -0.6108766748, 2.1930659258, 0.9091647433,
      0.3136661779, -0.6108766748, 1.4591897081, -1.1804364456, 0.3985316185,
      -1.0597154562, 2.1930659258, -1.1804364456, 3.3110327980, -1.1617171573,
      0.0083249547, 0.9091647433, 0.3985316185, -1.1617171573, 2.1174700695;

  Eigen::VectorXd b(5);
  b << -2.8634904630, -4.8216733374, -4.2958468309, -0.0877703331, -2.0223464006;

  double perturbationNorm = 1e-6;

  std::cout << "=== Analysis for matrix A1 ===\n";
  solveAndAnalyze(A1, b, perturbationNorm);

  std::cout << "=== Analysis for matrix A2 ===\n";
  solveAndAnalyze(A2, b, perturbationNorm);

  return 0;
}
