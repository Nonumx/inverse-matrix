#ifndef ELIMINATION_H
#define ELIMINATION_H
#include <cmath>
#include <stdexcept>
#include <vector>

#include "adjugate.h"

class EquationSystem {
  // EquationSystem solves only for the inverse
  // of a matrix, so initial matrix should be
  // (N, N) shape(N-order matrix).
 private:
  int row;
  std::vector<std::vector<double> > mat;

  int find_pivot(int);
  void interchange(const int&, const int&);
  void eliminate(int, int);

 public:
  EquationSystem(const std::vector<std::vector<double> >&, int);
  EquationSystem(const Matrix& m, int);
  ~EquationSystem() = default;
  void solve();
  double get_answer(int);
};

const double _eps = 1e-5;
double fgcd(double, double);
double flcm(double, double);

void example_elimination();
void solve_elimination(std::string filename);
#endif