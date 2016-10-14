#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace bnu = boost::numeric::ublas;
int solve_linear_system(const bnu::matrix<double>& A, bnu::vector<double>& x, bnu::vector<double>& y);
int invert_matrix(const bnu::matrix<double>& input, bnu::matrix<double>& inverse);

#endif // SOLVER_H_INCLUDED
