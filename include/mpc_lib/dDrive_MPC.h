#ifndef D_MPC
#define D_MPC

#include <Eigen-3.3/Eigen/QR>
#include <Eigen-3.3/Eigen/Core>
#include <vector>

struct MPCparams
{
    int N;
    double dt;

    double ref_cte;
    double ref_etheta;
    double ref_v;

    double MAX_OMEGA;
    double MAX_THROTTLE;
    double BOUND_VALUE;

    double W_CTE;
    double W_ETHETA;
    double W_VEL;
    double W_OMEGA;
    double W_ACC;
    double W_OMEGA_D;
    double W_ACC_D;

    MPCparams(void);
};

struct VarIndices
{
    size_t x_start, y_start, theta_start;
    size_t v_start, omega_start, acc_start;
    size_t cte_start, etheta_start;

    void setIndices(int N);
};

class MPC
{
public:
    MPC();

    virtual ~MPC();

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuatotions.
    std::vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, MPCparams params);
};

#endif