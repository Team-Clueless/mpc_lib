#include <mpc_lib/dDrive_MPC.h>
#include <Eigen-3.3/Eigen/QR>
#include <Eigen-3.3/Eigen/Core>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;

void VarIndices::setIndices(int N)
{
    x_start = 0;
    y_start = x_start + N;
    theta_start = y_start + N;
    v_start = theta_start + N;
    cte_start = v_start + N;
    etheta_start = cte_start + N;
    omega_start = etheta_start + N;
    acc_start = omega_start + N - 1;
}

MPCparams::MPCparams(void)
{
    BOUND_VALUE = 1.0e3;
}

struct MPCparams mpcParams;
struct VarIndices varIndices;

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
        for (int t = 0; t < mpcParams.N; t++)
        {
            fg[0] += mpcParams.W_CTE * CppAD::pow(vars[varIndices.cte_start + t] - mpcParams.ref_cte, 2);
            fg[0] += mpcParams.W_ETHETA * CppAD::pow(vars[varIndices.etheta_start + t] - mpcParams.ref_etheta, 2);
            fg[0] += mpcParams.W_VEL * CppAD::pow(vars[varIndices.v_start + t] - mpcParams.ref_v, 2);
        }
        for (int t = 0; t < mpcParams.N - 1; t++)
        {
            fg[0] += mpcParams.W_OMEGA * CppAD::pow(vars[varIndices.omega_start + t], 2);
            fg[0] += mpcParams.W_ACC * CppAD::pow(vars[varIndices.acc_start + t], 2);
        }
        // Smoother transitions (less jerks)
        for (int t = 0; t < mpcParams.N - 2; t++)
        {
            fg[0] += mpcParams.W_ACC_D * CppAD::pow(vars[varIndices.acc_start + t + 1] - vars[varIndices.acc_start + t], 2);
            fg[0] += mpcParams.W_OMEGA_D * CppAD::pow(vars[varIndices.omega_start + t + 1] - vars[varIndices.omega_start + t], 2);
        }
        //
        // Setup Constraints
        //

        // Initial constraints
        //
        // We add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`.
        // This bumps up the position of all the other values.
        fg[1 + varIndices.x_start] = vars[varIndices.x_start];
        fg[1 + varIndices.y_start] = vars[varIndices.y_start];
        fg[1 + varIndices.theta_start] = vars[varIndices.theta_start];
        fg[1 + varIndices.v_start] = vars[varIndices.v_start];
        fg[1 + varIndices.cte_start] = vars[varIndices.cte_start];
        fg[1 + varIndices.etheta_start] = vars[varIndices.etheta_start];

        // The rest of the constraints
        for (int t = 0; t < mpcParams.N - 1; t++)
        {

            // Time : T + 1
            AD<double> x1 = vars[varIndices.x_start + t + 1];
            AD<double> y1 = vars[varIndices.y_start + t + 1];
            AD<double> theta1 = vars[varIndices.theta_start + t + 1];
            AD<double> v1 = vars[varIndices.v_start + t + 1];
            AD<double> cte1 = vars[varIndices.cte_start + t + 1];
            AD<double> etheta1 = vars[varIndices.etheta_start + t + 1];

            // Time : T
            AD<double> x0 = vars[varIndices.x_start + t];
            AD<double> y0 = vars[varIndices.y_start + t];
            AD<double> theta0 = vars[varIndices.theta_start + t];
            AD<double> v0 = vars[varIndices.v_start + t];
            AD<double> cte0 = vars[varIndices.cte_start + t];
            AD<double> etheta0 = vars[varIndices.etheta_start + t];

            AD<double> w0 = vars[varIndices.omega_start + t];
            AD<double> a0 = vars[varIndices.acc_start + t];

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

            fg[2 + varIndices.x_start + t] = x1 - (x0 + v0 * CppAD::cos(theta0) * mpcParams.dt);
            fg[2 + varIndices.y_start + t] = y1 - (y0 + v0 * CppAD::sin(theta0) * mpcParams.dt);
            fg[2 + varIndices.theta_start + t] = theta1 - (theta0 + w0 * mpcParams.dt);
            fg[2 + varIndices.v_start + t] = v1 - (v0 + a0 * mpcParams.dt);

            fg[2 + varIndices.cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(etheta0) * mpcParams.dt));
            fg[2 + varIndices.etheta_start + t] = etheta1 - ((theta0 - traj_grad0) + w0 * mpcParams.dt);
        }
    }
};

MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, MPCparams params)
{
    mpcParams.N = params.N;
    mpcParams.dt = params.dt;

    mpcParams.ref_v = params.ref_v;
    mpcParams.ref_cte = params.ref_cte;
    mpcParams.ref_etheta = params.ref_etheta;

    mpcParams.MAX_OMEGA = params.MAX_OMEGA;
    mpcParams.MAX_THROTTLE = params.MAX_THROTTLE;

    mpcParams.W_CTE = params.W_CTE;
    mpcParams.W_ETHETA = params.W_ETHETA;
    mpcParams.W_VEL = params.W_VEL;
    mpcParams.W_OMEGA = params.W_OMEGA;
    mpcParams.W_ACC = params.W_ACC;
    mpcParams.W_OMEGA_D = params.W_OMEGA_D;
    mpcParams.W_ACC_D = params.W_ACC_D;

    varIndices.setIndices(mpcParams.N);

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
    size_t n_vars = 6 * mpcParams.N + 2 * (mpcParams.N - 1);
    // TODO: Set the number of constraints
    size_t n_constraints = 6 * mpcParams.N;

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++)
    {
        vars[i] = 0;
    }

    // set the initial variable values
    vars[varIndices.x_start] = x;
    vars[varIndices.y_start] = y;
    vars[varIndices.theta_start] = theta;
    vars[varIndices.v_start] = v;
    vars[varIndices.cte_start] = cte;
    vars[varIndices.etheta_start] = etheta;

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // TODO: Set lower and upper limits for variables.
    for (int i = 0; i < varIndices.omega_start; i++)
    {
        vars_lowerbound[i] = -mpcParams.BOUND_VALUE;
        vars_upperbound[i] = mpcParams.BOUND_VALUE;
    }
    for (int i = varIndices.omega_start; i < varIndices.acc_start; i++)
    {
        vars_lowerbound[i] = -mpcParams.MAX_OMEGA;
        vars_upperbound[i] = mpcParams.MAX_OMEGA;
    }
    for (int i = varIndices.acc_start; i < n_vars; i++)
    {
        vars_lowerbound[i] = -mpcParams.MAX_THROTTLE;
        vars_upperbound[i] = mpcParams.MAX_THROTTLE;
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
    constraints_lowerbound[varIndices.x_start] = x;
    constraints_lowerbound[varIndices.y_start] = y;
    constraints_lowerbound[varIndices.theta_start] = theta;
    constraints_lowerbound[varIndices.v_start] = v;
    constraints_lowerbound[varIndices.cte_start] = cte;
    constraints_lowerbound[varIndices.etheta_start] = etheta;

    constraints_upperbound[varIndices.x_start] = x;
    constraints_upperbound[varIndices.y_start] = y;
    constraints_upperbound[varIndices.theta_start] = theta;
    constraints_upperbound[varIndices.v_start] = v;
    constraints_upperbound[varIndices.cte_start] = cte;
    constraints_upperbound[varIndices.etheta_start] = etheta;

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

    result.push_back(solution.x[varIndices.omega_start]);
    result.push_back(solution.x[varIndices.acc_start]);

    // Add "future" solutions (where MPC is going)
    for (int i = 0; i < mpcParams.N - 1; ++i)
    {
        result.push_back(solution.x[varIndices.x_start + i + 1]);
        result.push_back(solution.x[varIndices.y_start + i + 1]);
    }

    return result;
}
