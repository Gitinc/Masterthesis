#include<cmath>
#include "MackeyGlass.h"
#include <fstream>
#include <stdexcept>
#include "iostream"
#include <iomanip>


mg_parameters::mg_parameters(double MGbeta, double MGgamma, int MGn, double MGtau){
    gamma = MGgamma;
    tau = MGtau;
    beta = MGbeta;
    n = MGn;
}


double deriv(double x, double  x_dot_tau,  mg_parameters MG)
{
    double dxdt = MG.beta * x_dot_tau / (1 + pow(x_dot_tau, MG.n))-(MG.gamma*x);
    return dxdt;
}
/*double mid_step(double older_state, double older_deriv, double newer_state, double newer_deriv){
    double ms_new = (1 + 2. * (0.5 - 0.0)/(1.0 - 0.0) ) * ((1.0 - 0.5)/(1.0 - 0.0))*((1.0 - 0.5)/(1.0 - 0.0)) * older_state
            +(0.5 - 0.0) * ((1.0 - 0.5)/(1.0 - 0.0))*((1.0 - 0.5)/(1.0 - 0.0)) * older_deriv
            +(1 + 2. * (0.5 - 1.0)/(0.0 - 1.0) ) * ((0.0 - 0.5)/(0.0 - 1.0)) * ((0.0 - 0.5)/(0.0 - 1.0)) * newer_state
            + (0.5 - 1.0) * ((0.0 - 0.5)/(0.0 - 1.0)) * ((0.0 - 0.5)/(0.0 - 1.0)) * newer_deriv;
    return ms_new;
}*/

//TODO: Add midsteps
double mackeyglass_rk4(double x_t, double x_t_minus_tau, double deltat, mg_parameters MG) {
    //double ms_new = mid_step(x_t_minus_tau);

    double k1 = deltat*deriv(x_t,          x_t_minus_tau, MG);
    double k2 = deltat*deriv(x_t+0.5*k1,   x_t_minus_tau, MG);
    double k3 = deltat*deriv(x_t+0.5*k2,   x_t_minus_tau, MG);
    double k4 = deltat*deriv(x_t+k3,       x_t_minus_tau, MG);
    double x_t_plus_deltat = (x_t + k1/6 + k2/3 + k3/3 + k4/6);
    return x_t_plus_deltat;
}

void mackey(double *X, double *T, int sample_n, mg_parameters MG, int interval) {
    double x0       = 1.2;		// initial condition: x(t=0)=x0
    double deltat   = 0.01;	    // time step size (which coincides with the integration step)


    double time = 0;
    int index = 1;
    int history_length = floor(MG.tau/deltat);
    double x_history[history_length];
    for (int i = 0; i < x_history[i]; ++i) x_history[i] = 0.0;
    double x_t = x0;
    double x_t_minus_tau, x_t_plus_deltat;

    // Set every value to the default value
    for (int i = 0; i < sample_n; i++) {
        if(i % interval == 0){
            X[(int) round(i/interval)] = x_t;
        }

        if (MG.tau == 0)
            x_t_minus_tau = 0.0;
        else
            x_t_minus_tau = x_history[index];


        x_t_plus_deltat = mackeyglass_rk4(x_t, x_t_minus_tau, deltat, MG);

        if (MG.tau != 0) {
            x_history[index] = x_t_plus_deltat;
            index = (index % history_length)+1;
        }

        time = time + deltat;
        if(i % interval == 0){
            T[(int) round(i/interval)] = time;
        }
        x_t = x_t_plus_deltat;
    }
}

void write_mackey_to_file(int sample_all, int interval){
    mg_parameters TestMG(0.2,0.1,10,17.0);
    double X[(int) round(sample_all/interval)];
    double T[(int) round(sample_all/interval)];
    for (int i = 0; i < (int) round(sample_all/interval); ++i) X[i] = 0.0;
    for (int i = 0; i < (int) round(sample_all/interval); ++i) T[i] = 0.0;
    mackey(X, T, sample_all, TestMG,interval);

    /*for (int i = 0; i < 15000; ++i) {
        std::cout << X[i*10] <<"\n";
    }*/
    std::ofstream myfile ("mackey.dat");
    if (myfile.is_open())
    {
        /*for(int i = 0; i <  (int) sample_all/interval; ++i){
            myfile << T[i*interval] << " ";
        }*/
        //myfile << "\n";
        for(int i = 0; i < (int) round(sample_all/interval); ++i){
            myfile << std::fixed << std::setprecision(15) << X[i] << "\n" ;
        }
        myfile.close();
    }
    else  throw std::invalid_argument( "File not open" );
}

/*int main(int argc, char*argv[]) {

    mg_parameters TestMG(0.2,0.1,10,17.0);
    int sample_all = 10000;    // total no. of samples, excluding the given initial condition
    double M[sample_all];
    double T[sample_all];
    for (int i = 0; i < sample_all; ++i) M[i] = 0.0;
    for (int i = 0; i < sample_all; ++i) T[i] = 0.0;

    mackey(M, T, sample_all, TestMG);
    for(int i = 0;i<sample_all;i++){
        std::cout << T[i] << ";" << M[i] << "\n";
    }
}*/
