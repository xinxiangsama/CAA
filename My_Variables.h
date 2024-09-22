// My_variables.h
#ifndef My_variables_H
#define My_variables_H
#include </home/xinxiangsama/lib/Eigen/Dense>
#include <vector>

double l = 100.0;
double Total_T = 50.0, x_min = -l,x_max = l, y_min = -l, y_max = l, Mx = 0.5, My = 0;
Eigen::MatrixXd Backward_coeff1(1, 7), Backward_coeff2(1, 7), Backward_coeff3(1, 7), Drp7(1, 7), b4(1, 4);
enum ComputeDirection {
    COMPUTE_UX,
    COMPUTE_UY,
    COMPUTE_BOTH
};

Eigen::MatrixXd fake(1,1);

std::vector<std::string> R_boundarys = {"left", "top", "bottom"};
std::vector<std::string> O_boundarys = {"right"};








#endif // My_variables_H