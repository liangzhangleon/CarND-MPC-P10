#include "MPC.h"
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

//Set the timestep length and duration
size_t N = 5;
double dt = 0.2;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 40 mph.
double ref_v = 70;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
 private:
  // The cost is stored in the first element of `fg`, any additions to the cost should be added to `fg[0]`.
  void SetCostF(ADvector& fg, const ADvector& vars) {
    fg[0] = 0;
    // The part of the cost based on the reference state.
    for (size_t t = 0; t < N; t++) {
      fg[0] += 10 * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    // Minimize the use of actuators.
    for (size_t t= 0; t < N - 1; t++) {
      fg[0] += 300 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }
    // Minimize the value gap between sequential actuations.
    for (size_t t= 0; t < N - 2; t++) {
      fg[0] += 5000 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
  }

  void SetConstraintsG(ADvector& fg, const ADvector& vars) {
    // Initial constraints
    // We add 1 to each of the starting indices due to cost being located at index 0 of `fg`.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (size_t t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(3.*coeffs[3]*x0 * x0 + 2.*coeffs[2]*x0 + coeffs[1]);

      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = y[t] - f(x[t])  + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      // The idea here is to move all Rhs of the equation to Lhs
      // Reformulate the model in the form as g(x) =0.
      fg[1 + x_start + t]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t]    = y1 - (y0 + v0 *  CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t]  = psi1 - (psi0 - v0 * delta0 / Lf * dt); //sign change here
      fg[1 + v_start + t]    = v1 - (v0 + a0 *  dt );
      fg[1 + cte_start + t]  = cte1 - ((y0 - f0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) - v0 *  delta0 / Lf * dt);
    }
  }

 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }
  // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
  void operator()(ADvector& fg, const ADvector& vars) {
    SetCostF(fg, vars);
    SetConstraintsG(fg, vars);
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

void MPC::InitVars(const Eigen::VectorXd &state, Dvector &vars) {
  for (size_t i = 0; i < vars.size(); i++) {
    vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start]    = state[0];
  vars[y_start]    = state[1];
  vars[psi_start]  = state[2];
  vars[v_start]    = state[3];
  vars[cte_start]  = state[4];
  vars[epsi_start] = state[5];
}

void MPC::SetVarsBound(Dvector &vars_lowerbound, Dvector &vars_upperbound) {
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = a_start; i < vars_lowerbound.size(); i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
}

void MPC::SetConstraintsBound(const Eigen::VectorXd &state, Dvector &constraints_lowerbound, Dvector &constraints_upperbound) {
  for (size_t i = 0; i < constraints_lowerbound.size(); i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start]    = state[0];
  constraints_lowerbound[y_start]    = state[1];
  constraints_lowerbound[psi_start]  = state[2];
  constraints_lowerbound[v_start]    = state[3];
  constraints_lowerbound[cte_start]  = state[4];
  constraints_lowerbound[epsi_start] = state[5];

  constraints_upperbound[x_start]    = constraints_lowerbound[x_start];
  constraints_upperbound[y_start]    = constraints_lowerbound[y_start];
  constraints_upperbound[psi_start]  = constraints_lowerbound[psi_start];
  constraints_upperbound[v_start]    = constraints_lowerbound[v_start];
  constraints_upperbound[cte_start]  = constraints_lowerbound[cte_start];
  constraints_upperbound[epsi_start] = constraints_lowerbound[epsi_start];
}


/*********************************************************************
* Ipopt solver is designed to solve such a optimization problem
*      min f(x) for x in R^n
*
*      s.t.   g_L <= g(x) <= g_U
*             x_L <=  x   <= x_U
* where f(x): R^n --> R is the objective function,
* g(x): R^n --> R^m are the constraint functions.
* g_L and g_U denote the lower and upper bounds on the constraints,
* x_L and x_U are the bounds on the variables x.
*
* To use Ipopt solver, one needs to prepare the following things
* 1. initialize x
* 2. set the lower and upper bound for x
* 3. set the lower and upper bound for g(x)
* 4. Define f(x) and g(x)
* ********************************************************************/
vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  // size_t i;

  // Set the number of model variables (includes both states and actuators).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // Set the number of constraints
  size_t n_constraints = 6 * N;

  // Declare and initialize value of the independent variables.
  // Should be 0 besides initial state.
  Dvector vars(n_vars);
  InitVars(state, vars);

  // Declare and set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  SetVarsBound(vars_lowerbound, vars_upperbound);

  // Declare and set Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  SetConstraintsBound(state, constraints_lowerbound, constraints_upperbound);


  // object that computes objective(also called as cost) and constraints
  FG_eval fg_eval(coeffs);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // Set a maximum time limit of the solver in seconds.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  vector<double> solved;
  solved.push_back(solution.x[delta_start]);
  solved.push_back(solution.x[a_start]);
  for (size_t i = 0; i < N; ++i) {
    solved.push_back(solution.x[x_start + i]);
    solved.push_back(solution.x[y_start + i]);
  }

  return solved;
}
