#ifndef DRP_H
#define DRP_H
#include "My_Variables.h"

#include </home/xinxiangsama/lib/Eigen/Dense>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath> 
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class DRP_solution
{   private:
    double dx, dy, dt;
    int nx, ny, nt ;
    double x, y;

    Eigen::MatrixXd U, V, P, rho; //velocity in x y , pressure, density
    Eigen::MatrixXd U_new, V_new, P_new, rho_new; //new velocity in x y , pressure, density
    Eigen::MatrixXd Udx, Vdx, Pdx, rhodx; //x_derivative of velocity in x y , pressure, density
    Eigen::MatrixXd Udy, Vdy, Pdy, rhody; //y_derivative of velocity in x y , pressure, density
    Eigen::MatrixXd Sin_theta, Cos_theta, V_theta, r;
    Eigen::MatrixXd Position_x, Position_y;
    Eigen::MatrixXd U_left_boundary, U_right_boundary, U_top_boundary, U_bottom_boundary;
    Eigen::MatrixXd V_left_boundary, V_right_boundary, V_top_boundary, V_bottom_boundary;
    Eigen::MatrixXd rho_left_boundary, rho_right_boundary, rho_top_boundary, rho_bottom_boundary;
    Eigen::MatrixXd P_left_boundary, P_right_boundary, P_top_boundary, P_bottom_boundary;
    Eigen::MatrixXd K1_U , K2_U , K3_U , K4_U ;
    Eigen::MatrixXd K1_V , K2_V , K3_V , K4_V ;
    Eigen::MatrixXd K1_P , K2_P , K3_P , K4_P ;
    Eigen::MatrixXd K1_rho , K2_rho , K3_rho , K4_rho ;
    


    public:
    DRP_solution(double dx_, double dy_, double dt_):dx(dx_), dy(dy_) , dt(dt_)
    {  
        nx = 2*l/dx + 1;
        ny = 2*l/dy + 1;
        nt = Total_T/dt + 1;
        U = V = P = rho = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        U_new = V_new = P_new = rho_new = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        Udx = Vdx = Pdx = rhodx = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        Udy = Vdy = Pdy = rhody = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        Position_x = Position_y = Eigen::MatrixXd::Zero(nx + 6, ny + 6);

        K1_U = K2_U = K3_U = K4_U = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        K1_V = K2_V = K3_V = K4_V = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        K1_P = K2_P = K3_P = K4_P = Eigen::MatrixXd::Zero(nx + 6, ny + 6);
        K1_rho = K2_rho = K3_rho = K4_rho = Eigen::MatrixXd::Zero(nx + 6, ny + 6);


        Backward_coeff1(0) = -2.192280339;Backward_coeff1(1) = 4.748611401;Backward_coeff1(2) = -5.108851915;Backward_coeff1(3) = 4.461567104;Backward_coeff1(4) = -2.833498741;Backward_coeff1(5) = 1.128328861;Backward_coeff1(6) = -0.203876371;
        Backward_coeff2(0) =  -0.209337622;Backward_coeff2(1) = -1.084875676;Backward_coeff2(2) = 2.147776050;Backward_coeff2(3) = -1.388928322;Backward_coeff2(4) = 0.768949766;Backward_coeff2(5) = -0.281814650;Backward_coeff2(6) = 0.048230454;
        Backward_coeff3(0) = 0.049041958;
        Backward_coeff3(1) = -0.468840357;
        Backward_coeff3(2) = -0.474760914;
        Backward_coeff3(3) =  1.273274737;
        Backward_coeff3(4) = -0.518484526;
        Backward_coeff3(5) = 0.166138533;
        Backward_coeff3(6) =  -0.026369431;
        Drp7(0) = -0.02084314277031176;
        Drp7(1) = 0.166705904414580469;
        Drp7(2) = -0.770882380518225552;
        Drp7(3) = 0.0;
        Drp7(4) = 0.770882380518225552;
        Drp7(5) = -0.166705904414580469;
        Drp7(6) = 0.02084314277031176;
        b4(0) = 2.3025580888;
        b4(1) = -2.4910075998;
        b4(2) = 1.5743409332;
        b4(3) = -0.3858914222;
        //initialize
        for(int i = 0; i < nx + 6; i++)
        {   x = (i - 3) * dx + x_min;

            for(int j = 0; j < ny + 6; j++)
            {
                y = (j - 3) * dy + y_min;
                Position_x(i, j) = x;
                Position_y(i, j) = y;
                U(i,j) = 0.04 * y * exp(-log(2.0) * (pow((x - 67.0), 2) + pow(y, 2)) / 25.0);
                V(i,j) = -0.04 * (x - 67) * exp(-log(2.0) * (pow((x - 67.0), 2) + pow(y, 2)) / 25.0);
                P(i,j) = exp(-log(2.0) * (pow(x, 2) + pow(y, 2)) / 9) ;
                rho(i,j) = exp(-log(2.0) * (pow(x, 2) + pow(y, 2)) / 9) + 0.1 *exp(-log(2.0) * (pow((x - 67.0), 2) + pow(y, 2)) / 25.0);
            }
        }
    }
    //边界及内部采用7点drp差分，边界外扩充的ghost points采用对应的向后差分
    void DRP_scheme(Eigen::MatrixXd& u_in, Eigen::MatrixXd& ux_out, Eigen::MatrixXd& uy_out, ComputeDirection direction)
    {  
        bool computeUx = (direction == COMPUTE_UX) || (direction == COMPUTE_BOTH);
        bool computeUy = (direction == COMPUTE_UY) || (direction == COMPUTE_BOTH);
        if(computeUx){
            //left bc======================================================================
            for(int n = 0; n < ny + 6; n ++)
            {
                for (int j = 0; j < 7; j ++)
                {
                    ux_out(0, n) += (Backward_coeff1(j) * u_in(0 + j, n))/dx ;
                }

                for (int j = -1; j < 6; j ++)
                {
                    ux_out(1, n) += (Backward_coeff2(j + 1) * u_in(1 + j, n))/dx ;
                }

                for (int j = -2; j < 5; j ++)
                {
                    ux_out(2, n) += (Backward_coeff3(j + 2) * u_in(2 + j, n))/dx ;
                }
            //right bc================================================================
                for (int j = 0; j < 7; j ++)
                {
                    ux_out(nx + 5, n) += ((-Backward_coeff1(j)) * u_in(nx + 5 - j, n))/dx ;
                }

                for (int j = -1; j < 6; j ++)
                {
                    ux_out(nx + 4, n) +=((-Backward_coeff2(j + 1)) * u_in(nx + 4 - j,n))/dx;
                }

                for (int j = -2; j < 5; j ++)
                {
                    ux_out(nx + 3, n) +=((-Backward_coeff2(j + 2)) * u_in(nx + 3 - j, n))/dx;
                }
            //interior points===========================================================
                for(int i = 3; i < nx +3; i++)
                {
                    for (int j = -3; j <= 3; j ++) 
                    {
                        ux_out(i, n) +=(1/dx) * (Drp7(j + 3) * u_in(i + j));
                    }   
                }
            }
        }

    if(computeUy){
        for(int i = 0; i < nx + 6; i ++)
        {
            //bottom bc======================================================================
            for (int k = 0; k < 7; k ++)
            {
                uy_out(i, 0) += (Backward_coeff1(k) * u_in(i, 0 + k))/dy ;
            }

            for (int k = -1; k < 6; k ++)
            {
                uy_out(i, 1) += (Backward_coeff2(k + 1) * u_in(i, 1 + k))/dy ;
            }

            for (int k = -2; k < 5; k ++)
            {
                uy_out(i, 2) += (Backward_coeff3(k + 2) * u_in(i, 2 + k))/dy ;
            }
            //top bc=====================================================================
            for (int k = 0; k < 7; k ++)
            {
                uy_out(i, ny + 5) += ((-Backward_coeff1(k)) * u_in(i, ny + 5 - k))/dy ;
            }

            for (int k = -1; k < 6; k ++)
            {
                uy_out(i, ny + 4) +=((-Backward_coeff2(k + 1)) * u_in(i, ny + 4 - k))/dy;
            }

            for (int k = -2; k < 5; k ++)
            {
                uy_out(i, ny + 3) +=((-Backward_coeff2(k + 2)) * u_in(i, ny + 3 - k))/dy;
            }
            //interior points===========================================================
            for(int j = 3; j < ny +3; j++)
            {
                for (int k = -3; k <= 3; k ++) 
                {
                    uy_out(i, j) +=(1/dy) * (Drp7(k + 3) * u_in(i, j + k));
                }   
            }
        }
    }
    }

    void domain_solve(Eigen::MatrixXd& u_out, const std::string& name)
    {   
        Eigen::MatrixXd E(nx + 6, ny + 6), F(nx + 6, ny + 6), Ex(nx + 6, ny + 6), Fy(nx + 6, ny + 6);

        // if(name == "rho"){E = Mx * rho + U; F = V;}
        // else if(name == "u"){E = Mx * rho + P; F = Eigen::MatrixXd::Zero(nx+6,nx+6);}
        // else if(name == "v"){E = Mx * V ; F = P;}
        // else if(name == "p"){E = Mx * P + U; F = V;}
        // else {cout<<"sucker! you submit meaningless name"<<endl;}
        if(name == "rho"){Ex = Mx * rhodx + Udx; Fy = Vdy;save_as_csv(Ex,"Ex_for_rho.csv");save_as_csv(Fy,"Fy_for_rho.csv");}
        else if(name == "u"){Ex = Mx * rhodx + Pdx; Fy = Eigen::MatrixXd::Zero(nx+6,nx+6);}
        else if(name == "v"){Ex = Mx * Vdx ; Fy = Pdy;}
        else if(name == "p"){Ex = Mx * Pdx + Udx; Fy = Vdy;}
        else {cout<<"sucker! you submit meaningless name"<<endl;} 

        
        // DRP_scheme(E, Ex, fake, COMPUTE_UX);
        // DRP_scheme(F, fake, Fy, COMPUTE_UY);
        for(int i = 3; i < nx + 3; i++)
        {
            for(int j = 3; j < ny + 3; j ++)
            {   
                u_out(i,j) = Ex(i, j) + Fy(i, j);
            }
        }
    }
   
   void Radiation_boundary_solve(const string& BC_name, const Eigen::MatrixXd& u, const Eigen::MatrixXd& partial_x, const Eigen::MatrixXd& partial_y, Eigen::MatrixXd K)
   {

    if(BC_name == "left")
    {   
        Sin_theta = Cos_theta = V_theta = r = Eigen::MatrixXd::Zero(1, ny + 6);
        for(int Dense = 1; Dense <=3 ; Dense ++){
            Get_position("left", Dense, Sin_theta, Cos_theta, V_theta, r);
            for(int i = 0; i < ny + 6; i ++)
            {
                K(Dense - 1, i) = V_theta(i) * (Cos_theta(i) * partial_x(Dense - 1,i) + Sin_theta(i) * partial_y(Dense - 1,i) + u(Dense - 1,i)/(2*r(i)));
            }
        }
    }

    else if(BC_name=="top")
    {
        Sin_theta = Cos_theta = V_theta = r = Eigen::MatrixXd::Zero(1, nx + 6);
        for(int Dense = 1; Dense <=3 ; Dense ++){
            Get_position("top", Dense, Sin_theta, Cos_theta, V_theta, r);
            for(int i = 0; i < nx + 6; i ++)
            {
                K(i, nx + 6 - Dense) = V_theta(i) * (Cos_theta(i) * partial_x(i,ny+6-Dense) + Sin_theta(i) * partial_y(i,ny+6-Dense) + u(i,ny+6-Dense)/(2*r(i)));
            }
        }
    }

    else if(BC_name == "bottom")
    {
        Sin_theta = Cos_theta = V_theta = r = Eigen::MatrixXd::Zero(1, nx + 6);
        for(int Dense = 1; Dense <=3 ; Dense ++){
            Get_position("bottom", Dense, Sin_theta, Cos_theta, V_theta, r);
            for(int i = 0; i < nx + 6; i ++)
            {
                K(i, Dense - 1) = V_theta(i) * (Cos_theta(i) * partial_x(i,Dense - 1) + Sin_theta(i) * partial_y(i,Dense - 1) + u(i,Dense - 1) /(2*r(i)));
            }
        }
    }
    // //右边界先试用辐射边界（本应该是出流边界，至少论文中是这样的）
    // else if(BC_name == "right")
    // {
    //     Sin_theta = Cos_theta = V_theta = r = Eigen::MatrixXd::Zero(1, nx + 6);
    //     for(int Dense = 1; Dense <=3 ; Dense ++){
    //         Get_position("right", Dense, Sin_theta, Cos_theta, V_theta, r);
    //         for(int i = 0; i < ny + 6; i ++)
    //         {
    //             K(nx + 6 - Dense, i) = V_theta(i) * (Cos_theta(i) * partial_x(nx+6-Dense,i) + Sin_theta(i) * partial_y(nx+6-Dense,i)  + u(nx+6-Dense,i) /(2*r(i)));
    //         }
    //     }
    // }

   }

    //out flow bc "right"
   void Outflow_boundary_solve(const string& BC_name, const string& name , Eigen::MatrixXd K)
   {
        if(BC_name == "right")
        {   
            Sin_theta = Cos_theta = V_theta = r = Eigen::MatrixXd::Zero(1, nx + 6);
            for(int Dense = 1; Dense <=3 ; Dense ++){
                Get_position("right", Dense, Sin_theta, Cos_theta, V_theta, r);
                for(int i = 0; i < ny + 6; i ++)
                {
                    if(name == "u")
                    {
                        K(nx + 6 - Dense, i) = Pdx(nx + 6 - Dense, i) + Mx * Udx(nx + 6 - Dense, i);
                    }
                    else if (name == "v")
                    {
                        K(nx + 6 - Dense, i) = Pdy(nx + 6 - Dense, i) + Mx * Vdx(nx + 6 - Dense, i);
                    }
                    else if (name == "p")
                    {
                        K(nx + 6 - Dense, i) = V_theta(i) * Cos_theta(i) * Pdx(nx + 6 - Dense, i) + V_theta(i) * Sin_theta(i) * Pdy(nx + 6 - Dense, i) + V_theta(i) * P(nx + 6 - Dense, i) /(2*r(i));
                    }

                    else if (name == "rho")
                    {
                        K(nx + 6 - Dense, i) = V_theta(i) * Cos_theta(i) * Pdx(nx + 6 - Dense, i) + V_theta(i) * Sin_theta(i) * Pdy(nx + 6 - Dense, i) + V_theta(i) * P(nx + 6 - Dense, i) /(2*r(i)) + Mx * (rhodx(nx + 6 - Dense, i) - Pdx( nx + 6 - Dense ,i));
                    }
                }
        }
        }

   }



    void Runge_kutta( const int& step)
    {   
        Eigen::MatrixXd U_temp(nx + 6, ny + 6), V_temp(nx + 6, ny + 6), P_temp(nx + 6, ny + 6), rho_temp(nx + 6, ny + 6);

        U_temp = U;
        V_temp = V;
        P_temp = P;
        rho_temp = rho;
        DRP_scheme(U_temp, Udx, Udy, COMPUTE_BOTH);
        DRP_scheme(V_temp, Vdx, Vdy, COMPUTE_BOTH);
        DRP_scheme(P_temp, Pdx, Pdy, COMPUTE_BOTH);
        DRP_scheme(rho_temp, rhodx, rhody, COMPUTE_BOTH);
        //first step y
        domain_solve(K1_U, "u");
        domain_solve(K1_V, "v");
        domain_solve(K1_P, "p");
        domain_solve(K1_rho, "rho");

        bc_condition(K1_U, K1_V, K1_P, K1_rho);

        //second step y + 1/2 * K1 * dt
        U_temp = U - 1/2 * K1_U * dt;
        V_temp = V - 1/2 * K1_V * dt;
        P_temp = P - 1/2 * K1_P * dt;
        rho_temp = rho - 1/2 * K1_rho * dt;
        DRP_scheme(U_temp, Udx, Udy, COMPUTE_BOTH);
        DRP_scheme(V_temp, Vdx, Vdy, COMPUTE_BOTH);
        DRP_scheme(P_temp, Pdx, Pdy, COMPUTE_BOTH);
        DRP_scheme(rho_temp, rhodx, rhody, COMPUTE_BOTH);
        domain_solve(K2_U, "u");
        domain_solve(K2_V, "v");
        domain_solve(K2_P, "p");
        domain_solve(K2_rho, "rho");

        bc_condition(K2_U, K2_V, K2_P, K2_rho);
        //third step y + 1/2 * K2 * dt
        U_temp = U - 1/2 * K2_U * dt;
        V_temp = V - 1/2 * K2_V * dt;
        P_temp = P - 1/2 * K2_P * dt;
        rho_temp = rho - 1/2 * K2_rho * dt;
        DRP_scheme(U_temp, Udx, Udy, COMPUTE_BOTH);
        DRP_scheme(V_temp, Vdx, Vdy, COMPUTE_BOTH);
        DRP_scheme(P_temp, Pdx, Pdy, COMPUTE_BOTH);
        DRP_scheme(rho_temp, rhodx, rhody, COMPUTE_BOTH);
        domain_solve(K3_U, "u");
        domain_solve(K3_V, "v");
        domain_solve(K3_P, "p");
        domain_solve(K3_rho, "rho");
        bc_condition(K3_U, K3_V, K3_P, K3_rho);
        //fourth step y + K3 * dt
        U_temp = U - K3_U * dt;
        V_temp = V - K3_V * dt;
        P_temp = P - K3_P * dt;
        rho_temp = rho - K3_rho * dt;
        DRP_scheme(U_temp, Udx, Udy, COMPUTE_BOTH);
        DRP_scheme(V_temp, Vdx, Vdy, COMPUTE_BOTH);
        DRP_scheme(P_temp, Pdx, Pdy, COMPUTE_BOTH);
        DRP_scheme(rho_temp, rhodx, rhody, COMPUTE_BOTH);
        domain_solve(K4_U, "u");
        domain_solve(K4_V, "v");
        domain_solve(K4_P, "p");
        domain_solve(K4_rho, "rho");
        bc_condition(K4_U, K4_V, K4_P, K4_rho);
        U_new = U - 1/6 * dt * (K1_U + 2 * K2_U + 2 * K3_U + K4_U);
        V_new = V - 1/6 * dt * (K1_V + 2 * K2_V + 2 * K3_V + K4_V);
        P_new = P - 1/6 * dt * (K1_P + 2 * K2_P + 2 * K3_P + K4_P);
        rho_new = rho - 1/6 * dt * (K1_rho + 2 * K2_rho + 2 * K3_rho + K4_rho);
    }

    void Adams_bathforth(const int& j)
    {
        Eigen::MatrixXd U_temp(nx + 6, ny + 6), V_temp(nx + 6, ny + 6), P_temp(nx + 6, ny + 6), rho_temp(nx + 6, ny + 6);
        U_temp = U;
        V_temp = V;
        P_temp = P;
        rho_temp = rho;
        DRP_scheme(U_temp, Udx, Udy, COMPUTE_BOTH);
        DRP_scheme(V_temp, Vdx, Vdy, COMPUTE_BOTH);
        DRP_scheme(P_temp, Pdx, Pdy, COMPUTE_BOTH);
        DRP_scheme(rho_temp, rhodx, rhody, COMPUTE_BOTH);
        //first step y
        domain_solve(K4_U, "u");
        domain_solve(K4_V, "v");
        domain_solve(K4_P, "p");
        domain_solve(K4_rho, "rho");

        bc_condition(K4_U, K4_V, K4_P, K4_rho);

        U_new = U - dt * (b4(0) * K4_U + b4(1) * K3_U + b4(2) * K2_U + b4(3) * K1_U);
        V_new = V - dt * (b4(0) * K4_V + b4(1) * K3_V + b4(2) * K2_V + b4(3) * K1_V);
        P_new = P - dt * (b4(0) * K4_P + b4(1) * K3_P + b4(2) * K2_P + b4(3) * K1_P);
        rho_new = rho - dt * (b4(0) * K4_rho + b4(1) * K3_rho + b4(2) * K2_rho + b4(3) * K1_rho);

        give_values(K1_U, K2_U, K3_U, K4_U);
        give_values(K1_V, K2_V, K3_V, K4_V);
        give_values(K1_P, K2_P, K3_P, K4_P);
        give_values(K1_rho, K2_rho, K3_rho, K4_rho);


    }

    void Runsim()
    {
        for(int step = 0; step < 500; step ++){
            // Adams_bathforth(step);
            Runge_kutta(step);
            U = U_new;
            V = V_new;
            P = P_new;
            rho = rho_new;

            cout<<step * dt<<endl;
        }
        save_as_csv(U,"U.csv");
        save_as_csv(V,"V.csv");
        save_as_csv(P,"P.csv");
        save_as_csv(rho,"rho.csv");

    }

    //Dense = 1,2,3
    void Get_position(const std::string& boundary, const int& Dense, Eigen::MatrixXd& sin, Eigen::MatrixXd& cos, Eigen::MatrixXd& v, Eigen::MatrixXd& r)
    {   
        //left boundary====================================================
        if(boundary == "left")
        {   
            for(int n = 0; n < ny + 6; n ++)
            {   
                r(n) = sqrt(pow(Position_x(Dense - 1, n), 2) + pow(Position_y(Dense - 1, n), 2));
                sin(n) = Position_y(Dense - 1, n)/r(n) ;
                cos(n) = Position_x(Dense - 1, n)/r(n);
                v(n) = Mx * cos(n) + sqrt(1 - pow((Mx * sin(n)), 2));
            }
        }
        //right boundary====================================================
        else if(boundary == "right")
        {
            for(int n = 0; n < ny + 6; n ++)
            {
                r(n) = sqrt(pow(Position_x(nx + 6 - Dense, n), 2) + pow(Position_y(nx + 6 - Dense, n), 2));
                sin(n) = Position_y(nx + 6 - Dense, n)/r(n) ;
                cos(n) = Position_x(nx + 6 - Dense, n)/r(n);
                v(n) = Mx * cos(n) + sqrt(1 - pow((Mx * sin(n)), 2));
            }
        }
        //top boundary====================================================
        else if(boundary == "top")
        {
            for(int n = 0; n < nx + 6; n ++)
            {
                r(n) = sqrt(pow(Position_x(n, ny + 6 - Dense), 2 )+ pow(Position_y(n, ny + 6 - Dense), 2));
                sin(n) = Position_y(n, ny + 6 - Dense)/r(n) ;
                cos(n) = Position_x(n, ny + 6 - Dense)/r(n);
                v(n) = Mx * cos(n) + sqrt(1 - pow((Mx * sin(n)), 2));
            }
        }
        //bottom boundary====================================================
        else if(boundary == "bottom")
        {
            for(int n = 0; n < nx + 6; n ++)
            {
                r(n) = sqrt(pow(Position_x(n, Dense - 1), 2) + pow(Position_y(n, Dense - 1), 2));
                sin(n) = Position_y(n, Dense - 1)/r(n) ;
                cos(n) = Position_x(n, Dense - 1)/r(n);
                v(n) = Mx * cos(n) + sqrt(1 - pow((Mx * sin(n)), 2));
            }
        }
        else
        {
            cout<<"sucker! you submit wrong BC_name!"<<endl;
        }
    }

    void save_as_csv(const Eigen::MatrixXd& matrix, const std::string& filename) const {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "无法打开文件: " << filename << std::endl;
            return;
        }
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j < matrix.cols() - 1) file << ","; // 在列之间添加逗号
            }
            if (i < matrix.rows() - 1) file << "\n"; // 在行之间添加换行符
        }

        file.close();
}

void give_values(Eigen::MatrixXd& m1, Eigen::MatrixXd& m2, Eigen::MatrixXd& m3, const Eigen::MatrixXd& m4 )
{
    m1 = m2;
    m2 = m3;
    m3 = m4;
}

void initialize_Kvalues(Eigen::MatrixXd& m1, Eigen::MatrixXd& m2, Eigen::MatrixXd& m3, const Eigen::MatrixXd& m4 )
{
    m1 = m4;
    m2 = m4;
    m3 = m4;
}

void bc_condition(Eigen::MatrixXd& K_U,Eigen::MatrixXd& K_V, Eigen::MatrixXd& K_P, Eigen::MatrixXd& K_rho )
{
    for (const string& bc : R_boundarys) 
    {
        Radiation_boundary_solve(bc, U, Udx, Udy, K_U);
        Radiation_boundary_solve(bc, V, Vdx, Vdy, K_V);
        Radiation_boundary_solve(bc, P, Pdx, Pdy, K_P);
        Radiation_boundary_solve(bc, rho, rhodx, rhody, K_rho);
    }    
        

    for (const string& bc : O_boundarys) 
    {
        Outflow_boundary_solve(bc, "u", K_U);
        Outflow_boundary_solve(bc, "v", K_V);
        Outflow_boundary_solve(bc, "p", K_P);
        Outflow_boundary_solve(bc, "rho", K_rho);
    }
}

    
};




#endif // DRP_H