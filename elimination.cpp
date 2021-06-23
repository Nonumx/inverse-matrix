#include "elimination.h"

int rank, proc;
int q, r;
int *counts; //每个进程传输的数据个数
int *displs; //偏移量

EquationSystem::EquationSystem(const std::vector<std::vector<double>> &x,
                               int current_column)
{
  this->row = x.size();
  for (auto &v : x)
  {
    if (v.size() != x.size())
    {
      throw std::logic_error("Matrix should be N-order.");
    }
  }

  for (int i = 0; i < x.size(); i++)
  {
    int equation_index = this->mat.size();
    this->mat.push_back(std::vector<double>());
    for (int k = 0; k < x.size(); k++)
    {
      mat[equation_index].push_back(x[i][k]);
    }
    if (i == current_column)
    {
      // diagonal
      mat[equation_index].push_back(1.0);
    }
    else
    {
      mat[equation_index].push_back(0.0);
    }
  }
}

EquationSystem::EquationSystem(const Matrix &m, int current_column)
{
  // Syntactic sugar
  *this = EquationSystem(m.mat, current_column);
}

int EquationSystem::find_pivot(int row_entry)
{
  int pivot_row = row_entry;
  int leftmost_nonzero = mat[0].size();
  volatile bool flag = false;
#pragma omp parallel for num_threads(8) private(mat, i, j) shared(flag) firstprivate(pivot_row, leftmost_nonzero)
  for (int i = row_entry; i < this->row && !flag; i++)
  {
    // elimination starts from row 0, and the submatrix to be
    // eliminated will be echelon form. (Numbers before `row_entry`
    // are all zero)
    for (int j = row_entry; j < mat[i].size(); j++)
    {
      if (flag)
        continue;
      if (fabs(mat[i][j]) > _eps)
      {
        if (leftmost_nonzero > j)
        {
          pivot_row = i;
          leftmost_nonzero = j;
        }
        flag = true;
      }
    }
  }
  return pivot_row;
}

void EquationSystem::interchange(const int &row_a, const int &row_b)
{
  std::swap(this->mat[row_a], this->mat[row_b]);
}

double fgcd(double x, double y)
{
  double r = x - floor(x / y) * y;
  while (fabs(r) > _eps)
  {
    x = y;
    y = r;
    r = x - floor(x / y) * y;
  }
  return y;
}

double flcm(double x, double y) { return x / fgcd(x, y) * y; }

void EquationSystem::eliminate(int pivot_row, int eliminate_row)
{
  // echelon form matrix, pivot_row has be
  // placed to the first row of submatrix
  double x = this->mat[pivot_row][pivot_row];
  double y = this->mat[eliminate_row][pivot_row];
  if (fabs(y) < _eps)
  {
    return;
  }
  double lcm = flcm(x, y);
  double factor_x = lcm / x;
  double factor_y = lcm / y;

  for (int i = pivot_row; i < this->mat[eliminate_row].size(); i++)
  {
    this->mat[eliminate_row][i] = this->mat[eliminate_row][i] * factor_y -
                                  this->mat[pivot_row][i] * factor_x;
  }
}

void EquationSystem::solve()
{
  for (int r = 0; r < this->row; r++)
  {
    int pivot_row = find_pivot(r);
    interchange(r, pivot_row);
    for (int i = r + 1; i < this->row; i++)
    {
      eliminate(r, i);
    }
  }
  for (int r = 0; r < this->row; r++)
  {
    double scale = this->mat[r][r];
    for (int i = r; i < this->mat[r].size(); i++)
    {
      this->mat[r][i] /= scale;
    }
  }

  for (int r = this->row - 1; r >= 0; r--)
  {
    double ans = *(this->mat[r].end() - 1);
    for (int i = r + 1; i < this->mat[r].size() - 1; i++)
    {
      ans -= *(this->mat[i].end() - 1) * this->mat[r][i];
      this->mat[r][i] = 0.0;
    }
    *(this->mat[r].end() - 1) = ans;
  }
}

double EquationSystem::get_answer(int i) { return *(mat[i].end() - 1); }

void example_elimination()
{
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
  std::cout << m;
  std::cout << "\n";

  Matrix inv_m(3, 3);

  for (int i = 0; i < m.col; i++)
  {
    EquationSystem eq(m, i);
    eq.solve();
    for (int j = 0; j < m.row; j++)
    {
      inv_m(j, i) = eq.get_answer(j);
    }
  }

  std::cout << "inverse:\n";
  std::cout << inv_m;
  std::cout << "\n";

  m = matmul(m, inv_m);
  std::cout << "final:\n";
  std::cout << m;
  std::cout << "\n";
}

void solve_elimination(std::string filename)
{
  // "/home/nonumx/CLionProjects/inv/matrix/arc130/arc130.mtx"
  double starttime, endtime;
  int row, column, nz;
  int i, j;
  std::ifstream ifs;
  ifs.open(filename);
  if (!rank)
  {
    if (!ifs.is_open())
    {
      std::cout << "failed to open file: " << filename << "\n";
      return;
    }
    ifs >> row >> column >> nz;
  }

  MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&column, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<double> intermediate(row * column, 0);

  Matrix m(row, column);

  if (!rank)
  {
    double v;
    while (ifs >> i >> j >> v)
    {
      i--;
      j--;
      m(i, j) = v;
      intermediate[i * column + j] = v;
    }
  }
  
  MPI_Bcast(&intermediate[0], row * column, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (rank)
  {
    for (i = 0; i < row; ++i)
      for (j = 0; j < column; ++j)
        m(i, j) = intermediate[i * column + j];
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);

  counts = (int *)malloc(proc * sizeof(int));
  displs = (int *)malloc(proc * sizeof(int));
  q = row / proc; //计算偏移量
  r = row % proc;
  for (i = 0; i < proc; ++i)
  {
    counts[i] = (q + ((i < r) ? 1 : 0)) * row;
    displs[i] = (q * i + ((i < r) ? i : r)) * row;
  }

  if (!rank)
    starttime = omp_get_wtime();

  Matrix inv_m(row, column);
  std::vector<double> intermediate_(counts[rank], 0);
  std::vector<double> result(row * column, 0);

  for (int i = 0; i < counts[rank] / row; i++)
  {
    EquationSystem eq(m, i + displs[rank] / row);
    eq.solve();
    for (int j = 0; j < m.row; j++)
    {
      intermediate_[i * column + j] = eq.get_answer(j);
    }
  }

  MPI_Gatherv(&intermediate_[0], counts[rank], MPI_DOUBLE, &result[0], counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (!rank)
  {
    for (i = 0; i < row; ++i)
      for (j = 0; j < column; ++j)
        inv_m(j, i) = result[i * column + j];
    endtime = omp_get_wtime();
  }

  if (!rank)
  {
    std::cout << "Use time: " << endtime - starttime << std::endl;
    m = matmul(m, inv_m);
    std::cout << "check:\n";
    bool flag = true;
    for (int i = 0; i < m.row; i++)
    {
      for (int j = 0; j < m.col; j++)
      {
        if (i == j)
        {
          if (fabs(m(i, j) - 1.0) > _eps)
          {
            flag = false;
            break;
          }
        }
        else
        {
          if (fabs(m(i, j)) > _eps)
          {
            flag = false;
            break;
          }
        }
      }
      if (!flag)
        break;
    }
    if (flag)
    {
      std::cout << "correct\n";
    }
    else
    {
      std::cout << "error\n";
    }
  }
}