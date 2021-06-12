#include "adjugate.h"

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

Matrix cofactor(const Matrix &m, int select_row, int select_col) {
  Matrix res(m.row - 1, m.col - 1);
  for (int i = 0; i < m.row; i++) {
    for (int j = 0; j < m.col; j++) {
      if (i == select_row || j == select_col) {
        continue;
      }

      if (i >= select_row && j < select_col) {
        res(i - 1, j) = m(i, j);
      } else if (i < select_row && j >= select_col) {
        res(i, j - 1) = m(i, j);
      } else if (i >= select_row && j >= select_col) {
        res(i - 1, j - 1) = m(i, j);
      } else {
        res(i, j) = m(i, j);
      }
    }
  }
  return res;
}
Matrix T(const Matrix &m) {
  Matrix t(m.col, m.row);
  for (int i = 0; i < m.row; i++) {
    for (int j = 0; j < m.col; j++) {
      t(j, i) = m(i, j);
    }
  }
  return t;
}

double det(const Matrix &m) {
  double ans = 0.0;

  std::vector<int> idx(m.row);
  for (int i = 0; i < m.row; i++) {
    idx[i] = i;
  }

  do {
    int sgn = 0;
    for (int i = 0; i < m.row; i++) {
      for (int j = i + 1; j < m.row; j++) {
        if (idx[j] < idx[i]) sgn++;
      }
    }

    if (sgn & 1) {
      sgn = -1;
    } else {
      sgn = 1;
    }

    double prod = 1.0;
    for (int i = 0; i < m.row; i++) {
      prod *= m(i, idx[i]);
    }
    prod *= sgn;

    ans += prod;
  } while (std::next_permutation(idx.begin(), idx.end()));
  return ans;
}

Matrix adjugate(const Matrix &m) {
  Matrix ans(m.row, m.col);
  for (int i = 0; i < m.row; i++) {
    for (int j = 0; j < m.col; j++) {
      ans(i, j) = det(cofactor(m, i, j));
      if ((i + j) & 1) {
        ans(i, j) = -ans(i, j);
      }
    }
  }
  ans = T(ans);
  return ans;
}

Matrix divide(const Matrix &m, double v) {
  Matrix ans(m.row, m.col);
  for (int i = 0; i < m.row; i++) {
    for (int j = 0; j < m.col; j++) {
      ans(i, j) = m(i, j) / v;
    }
  }
  return ans;
}

Matrix inverse(const Matrix &m) {
  Matrix adj = adjugate(m);
  double det_adj = det(m);
  Matrix inv_m = divide(adj, det_adj);
  return inv_m;
}

void example1() {
  Matrix m(3, 3);
  m(0, 0) = 2;
  m(0, 1) = 3;
  m(0, 2) = 1;
  m(1, 0) = 3;
  m(1, 1) = 4;
  m(1, 2) = 1;
  m(2, 0) = 3;
  m(2, 1) = 7;
  m(2, 2) = 2;
  std::cout << "Matrix m is:\n";
  std::cout << m << std::endl;

  Matrix inv_m = inverse(m);
  std::cout << "Inverse of m:\n";
  std::cout << inv_m << "\n";

  m = matmul(m, inv_m);
  std::cout << "Check:\n";
  std::cout << m << "\n";
}

void solve(std::string filename) {
  // "/home/nonumx/CLionProjects/inv/matrix/arc130/arc130.mtx"
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    return;
  }

  int row, column, nz;
  ifs >> row >> column >> nz;
  Matrix m(row, column);

  int i, j;
  double v;
  while (ifs >> i >> j >> v) {
    i--;
    j--;
    m(i, j) = v;
  }

  Matrix inv_m = inverse(m);

  m = matmul(m, inv_m);
  std::cout << "check:\n";
  std::cout << m;
}