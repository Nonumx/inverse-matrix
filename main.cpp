#include "elimination.h"

int main() {
  MPI_Init(nullptr, nullptr);
  solve_elimination("../data/nasa1824/nasa1824.mtx");
  MPI_Finalize();
  return 0;
}