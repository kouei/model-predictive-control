#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using Eigen::VectorXd;
using namespace std;

// TODO: Set the timestep length and duration
auto N = 12;
auto dt = 0.11;
auto latency = 0.1;

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
const auto Lf = 2.67;

auto ref_v = 100.0;
auto x_start = 0;
auto y_start = x_start + N;
auto psi_start = y_start + N;
auto v_start = psi_start + N;
auto cte_start = v_start + N;
auto epsi_start = cte_start + N;
auto delta_start = epsi_start + N;
auto a_start = delta_start + N - 1;

class FG_eval
{
public:
  // Fitted polynomial coefficients
  VectorXd coeffs;
  FG_eval(VectorXd coeffs_) : coeffs(coeffs_) { }

  using ADDouble = AD<double>;
  using ADvector = CPPAD_TESTVECTOR(ADDouble);

  template<typename T>
  T square(const T & a) const
  {
    return a * a;
  }

  template<typename T>
  T cube(const T & a) const
  {
    return a * a * a;
  }

  vector<double> normalize_weight(const vector<double> & weight)
  {
    vector<double> result;

    auto sum = 0.0;
    for(auto w : weight) sum += w;
    for(auto w : weight) result.push_back(w / sum);

    return result;
  }

  void operator()(ADvector & fg, const ADvector & vars)
  {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // normalized weight for each cost
    static auto weight = normalize_weight(vector<double> { 4000, 4000, 0.6, 10, 10, 500, 300, 50 });

    // Reference State Cost
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.

    // minimize the Cross Track Error and Steer Angle Error
    // minimize the gap between current velocity and reference velocity
    for (auto i = 0; i < N; ++i)
    {
      fg[0] += weight[0] * square(vars[cte_start + i]);
      fg[0] += weight[1] * square(vars[epsi_start + i]);
      fg[0] += weight[2] * square(vars[v_start + i] - ref_v);
    }

    // minimize usage of actuators
    for (auto i = 0; i < N - 1; ++i)
    {
      fg[0] += weight[3] * square(vars[delta_start + i]);
      fg[0] += weight[4] * square(vars[a_start + i]);
      
      // penalty for speed times steer
      // this will make sure the vehicle won't turn around in a high speed
      fg[0] += weight[5] * square(vars[delta_start + i] * vars[v_start + i]);
    }

    // smooth the trajectory, minimize the gap between adjacent usage of actuators
    for (auto i = 0; i < N - 2; i++)
    {
      fg[0] += weight[6] * square(vars[delta_start + i + 1] - vars[delta_start + i]);
      fg[0] += weight[7] * square(vars[a_start + i + 1] - vars[a_start + i]);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (auto t = 1; t < N; ++t)
    {
      auto x1 = ADDouble(vars[x_start + t]);
      auto x0 = ADDouble(vars[x_start + t - 1]);
      auto y1 = ADDouble(vars[y_start + t]);
      auto y0 = ADDouble(vars[y_start + t - 1]);
      auto psi1 = ADDouble(vars[psi_start + t]);
      auto psi0 = ADDouble(vars[psi_start + t - 1]);
      auto v1 = ADDouble(vars[v_start + t]);
      auto v0 = ADDouble(vars[v_start + t - 1]);
      auto cte1 = ADDouble(vars[cte_start + t]);
      auto cte0 = ADDouble(vars[cte_start + t - 1]);
      auto epsi1 = ADDouble(vars[epsi_start + t]);
      auto epsi0 = ADDouble(vars[epsi_start + t - 1]);
      auto a = ADDouble(vars[a_start + t - 1]);
      auto delta = ADDouble(vars[delta_start + t - 1]);
      if (t > 1)
      { // use previous actuations (to account for latency)
        a = vars[a_start + t - 2];
        delta = vars[delta_start + t - 2];
      }
      auto f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * square(x0) + coeffs[3] * cube(x0);
      auto psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * square(x0));

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.

      // TODO: Setup the rest of the model constraints
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 - v0 / Lf * delta * dt);
      fg[1 + v_start + t] = v1 - (v0 + a * (dt - latency));
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) - v0 / Lf * delta * (dt - latency));
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
  auto ok = true;
  using Dvector = CPPAD_TESTVECTOR(double);

  auto x = state[0];
  auto y = state[1];
  auto psi = state[2];
  auto v = state[3];
  auto cte = state[4];
  auto epsi = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  auto n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  auto n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  auto vars = Dvector(n_vars);
  for (auto i = 0; i < n_vars; i++)
  {
    vars[i] = 0;
  }

  auto vars_lowerbound = Dvector(n_vars);
  auto vars_upperbound = Dvector(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (auto i = 0; i < delta_start; ++i)
  {
    vars_lowerbound[i] = -numeric_limits<double>::infinity();
    vars_upperbound[i] = numeric_limits<double>::infinity();
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (auto i = delta_start; i < a_start; ++i)
  {
    vars_lowerbound[i] = -M_PI * 25 / 180;
    vars_upperbound[i] = M_PI * 25 / 180;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (auto i = a_start; i < n_vars; ++i)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (auto i = 0; i < n_constraints; ++i)
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  auto fg_eval = FG_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
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
  //std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (auto i = 0; i < N - 1; ++i)
  {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
