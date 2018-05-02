#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

using namespace std;
using CppAD::AD;

class MPC {
 public:
  typedef CPPAD_TESTVECTOR(double) Dvector;
 private:
  void InitVars(const Eigen::VectorXd &state, Dvector &vars);

  void SetVarsBound(Dvector &vars_lowerbound, Dvector &vars_upperbound);

  void SetConstraintsBound(const Eigen::VectorXd &state, Dvector &constraints_lowerbound, Dvector &constraints_upperbound);
 public:
  MPC();
  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
