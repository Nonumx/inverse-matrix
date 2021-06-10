#include <iostream>
#include <vector>

struct Matrix {
  int col;
  int row;
  std::vector<std::vector<double> > mat;

  Matrix(int row, int col) {
    this->row = row;
    this->col = col;

    for (int i = 0; i < row; i++) {
      mat.emplace_back();
      for (int j = 0; j < col; j++) {
        mat[i].push_back(0.0);
      }
    }
  }

  friend std::ostream &operator<<(std::ostream &os, Matrix &m) {
    for (int i = 0; i < m.row; i++) {
      for (int j = 0; j < m.col; j++) {
        os << m(i, j);
        if (j != m.col - 1) {
          os << " ";
        }
      }
      os << std::endl;
    }
    return os;
  }

  double operator()(int i, int j) const { return mat[i][j]; }

  double &operator()(int i, int j) { return mat[i][j]; }
};

Matrix matmul(const Matrix &a, const Matrix &b) {
  Matrix r(a.row, b.col);
  for (int i = 0; i < a.row; i++) {
    for (int j = 0; j < b.col; j++) {
      for (int k = 0; k < a.col; k++) {
        r(i, j) = r(i, j) + a(i, k) * b(k, j);
      }
    }
  }
  return r;
}