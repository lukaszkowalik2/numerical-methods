#ifndef HELPER_H
#define HELPER_H

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

void printMatrix(const std::vector<std::vector<double>>& matrix) {
  std::cout << "[" << std::endl;
  for (size_t i = 0; i < matrix.size(); i++) {
    std::cout << "  [";
    for (size_t j = 0; j < matrix[i].size(); j++) {
      std::cout << matrix[i][j];
      if (j < matrix[i].size() - 1) std::cout << ", ";
    }
    std::cout << "]";
    if (i < matrix.size() - 1) std::cout << ",";
    std::cout << std::endl;
  }
  std::cout << "]" << std::endl;
}

void printVector(const std::vector<double>& vector) {
  std::cout << "[";
  for (size_t i = 0; i < vector.size(); i++) {
    std::cout << vector[i];
    if (i < vector.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
}

Eigen::MatrixXd convertMatrixToEigen(const std::vector<std::vector<double>>& matrix) {
  Eigen::MatrixXd eigenMatrix(matrix.size(), matrix[0].size());
  for (size_t i = 0; i < matrix.size(); i++) {
    for (size_t j = 0; j < matrix[i].size(); j++) {
      eigenMatrix(i, j) = matrix[i][j];
    }
  }
  return eigenMatrix;
}

std::vector<std::vector<double>> convertEigenToMatrix(const Eigen::MatrixXd& eigenMatrix) {
  std::vector<std::vector<double>> matrix(eigenMatrix.rows(), std::vector<double>(eigenMatrix.cols()));
  for (int i = 0; i < eigenMatrix.rows(); i++) {
    for (int j = 0; j < eigenMatrix.cols(); j++) {
      matrix[i][j] = eigenMatrix(i, j);
    }
  }
  return matrix;
}

Eigen::VectorXd convertVectorToEigen(const std::vector<double>& vector) {
  Eigen::VectorXd eigenVector(vector.size());
  for (size_t i = 0; i < vector.size(); i++) {
    eigenVector(i) = vector[i];
  }
  return eigenVector;
}

std::vector<double> convertEigenToVector(const Eigen::VectorXd& eigenVector) {
  std::vector<double> vector(eigenVector.size());
  for (int i = 0; i < eigenVector.size(); i++) {
    vector[i] = eigenVector(i);
  }
  return vector;
}

void printEigenMatrix(const Eigen::MatrixXd& matrix) {
  std::cout << matrix << std::endl;
}

void printEigenVector(const Eigen::VectorXd& vector) {
  std::cout << vector << std::endl;
}

#endif  // HELPER_H
