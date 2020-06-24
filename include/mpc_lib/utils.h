#ifndef MPC_U
#define MPC_U

#include <cmath>
#include <Eigen-3.3/Eigen/QR>
#include <Eigen-3.3/Eigen/Core>

double pi(void);

double deg2rad(double x);
double rad2deg(double x);

double polyeval(Eigen::VectorXd coeffs, double x);
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);

#endif
