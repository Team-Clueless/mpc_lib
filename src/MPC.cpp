#include <mpc_lib/MPC.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <Eigen-3.3/Eigen/Core>

using CppAD::AD;

// Number of timesteps
int N = 12;
// Duration of each timestep
double dt = 0.06;
// Reference states
double ref_cte = 0;
double ref_etheta = 0;
double ref_v = 12;

// setup weights for cost function
const double W_CTE = 5000;
const double W_ETHETA = 3000;
const double W_VEL = 100;
const double W_OMEGA = 5;
const double W_ACC = 5;
const double W_OMEGA_D = 10;
const double W_ACC_D = 500;

// set bounds for actuator values
const double MAX_OMEGA = 3.0;
const double MAX_THROTTLE = 2.0;
const double BOUND_VALUE = 1.0e3; // every other

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our life easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t theta_start = y_start + N;
size_t v_start = theta_start + N;
size_t cte_start = v_start + N;
size_t etheta_start = cte_start + N;
size_t omega_start = etheta_start + N;
size_t acc_start = omega_start + N - 1;

class FG_eval
{
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;
    FG_eval(Eigen::VectorXd coeffs)
    {
        this->coeffs = coeffs;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector &fg, const ADvector &vars)
    {
        // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)

        // The cost is stored is the first element of `fg`.
        // Any additions to the cost should be added to `fg[0]`.
        fg[0] = 0;

        // Reference State Cost
        for (int t = 0; t < N; t++)
        {
            fg[0] += W_CTE * CppAD::pow(vars[cte_start + t] - ref_cte, 2);
            fg[0] += W_ETHETA * CppAD::pow(vars[etheta_start + t] - ref_etheta, 2);
            fg[0] += W_VEL * CppAD::pow(vars[v_start + t] - ref_v, 2);
        }
        for (int t = 0; t < N - 1; t++)
        {
            fg[0] += W_OMEGA * CppAD::pow(vars[omega_start + t], 2);
            fg[0] += W_ACC * CppAD::pow(vars[acc_start + t], 2);
        }
        // Smoother transitions (less jerks)
        for (int t = 0; t < N - 2; t++)
        {
            fg[0] += W_ACC_D * CppAD::pow(vars[acc_start + t + 1] - vars[acc_start + t], 2);
            fg[0] += W_OMEGA_D * CppAD::pow(vars[omega_start + t + 1] - vars[omega_start + t], 2);
        }
        //
        // Setup Constraints
        //

        // Initial constraints
        //
        // We add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`.
        // This bumps up the position of all the other values.
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + theta_start] = vars[theta_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + etheta_start] = vars[etheta_start];

        // The rest of the constraints
        for (int t = 0; t < N - 1; t++)
        {

            // Time : T + 1
            AD<double> x1 = vars[x_start + t + 1];
            AD<double> y1 = vars[y_start + t + 1];
            AD<double> theta1 = vars[theta_start + t + 1];
            AD<double> v1 = vars[v_start + t + 1];
            AD<double> cte1 = vars[cte_start + t + 1];
            AD<double> etheta1 = vars[etheta_start + t + 1];

            // Time : T
            AD<double> x0 = vars[x_start + t];
            AD<double> y0 = vars[y_start + t];
            AD<double> theta0 = vars[theta_start + t];
            AD<double> v0 = vars[v_start + t];
            AD<double> cte0 = vars[cte_start + t];
            AD<double> etheta0 = vars[etheta_start + t];

            AD<double> w0 = vars[omega_start + t];
            AD<double> a0 = vars[acc_start + t];

            AD<double> f0 = 0.0;
            for (int i = 0; i < coeffs.size(); i++)
            {
                f0 += coeffs[i] * CppAD::pow(x0, i);
            }

            AD<double> traj_grad0 = 0.0;
            for (int i = 1; i < coeffs.size(); i++)
            {
                traj_grad0 += i * coeffs[i] * CppAD::pow(x0, i - 1);
            }
            traj_grad0 = CppAD::atan(traj_grad0);

            // Here's `x` to get you started.
            // The idea here is to constraint this value to be 0.
            //
            // NOTE: The use of `AD<double>` and use of `CppAD`!
            // This is also CppAD can compute derivatives and pass
            // these to the solver.

            fg[2 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(theta0) * dt);
            fg[2 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(theta0) * dt);
            fg[2 + theta_start + t] = theta1 - (theta0 + w0 * dt);
            fg[2 + v_start + t] = v1 - (v0 + a0 * dt);

            fg[2 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(etheta0) * dt));
            fg[2 + etheta_start + t] = etheta1 - ((theta0 - traj_grad0) + w0 * dt);
        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    const double x = state[0];
    const double y = state[1];
    const double theta = state[2];
    const double v = state[3];
    const double cte = state[4];
    const double etheta = state[5];

    // TODO: Set the number of model variables (includes both states and inputs).
    // For example: If the state is a 4 element vector, the actuators is a 2
    // element vector and there are 10 timesteps. The number of variables is:
    //
    // 4 * 10 + 2 * 9
    size_t n_vars = 6 * N + 2 * (N - 1);
    // TODO: Set the number of constraints
    size_t n_constraints = 6 * N;

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++)
    {
        vars[i] = 0;
    }

    // set the initial variable values
    vars[x_start] = x;
    vars[y_start] = y;
    vars[theta_start] = theta;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[etheta_start] = etheta;

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // TODO: Set lower and upper limits for variables.
    for (int i = 0; i < omega_start; i++)
    {
        vars_lowerbound[i] = -BOUND_VALUE;
        vars_upperbound[i] = BOUND_VALUE;
    }
    for (int i = omega_start; i < acc_start; i++)
    {
        vars_lowerbound[i] = -MAX_OMEGA;
        vars_upperbound[i] = MAX_OMEGA;
    }
    for (int i = acc_start; i < n_vars; i++)
    {
        vars_lowerbound[i] = -MAX_THROTTLE;
        vars_upperbound[i] = MAX_THROTTLE;
    }
    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);

    for (int i = 0; i < n_constraints; i++)
    {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[theta_start] = theta;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[etheta_start] = etheta;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[theta_start] = theta;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[etheta_start] = etheta;

    // object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    //
    // NOTE: You don't have to worry about these options
    //
    // options for IPOPT solver
    std::string options;
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
    std::cout << "COST  : " << solution.obj_value << std::endl;

    // TODO: Return the first actuator values. The variables can be accessed with
    // `solution.x[i]`.
    std::vector<double> result;

    result.push_back(solution.x[omega_start]);
    result.push_back(solution.x[acc_start]);

    // Add "future" solutions (where MPC is going)
    for (int i = 0; i < N - 1; ++i)
    {
        result.push_back(solution.x[x_start + i + 1]);
        result.push_back(solution.x[y_start + i + 1]);
    }

    return result;
}
