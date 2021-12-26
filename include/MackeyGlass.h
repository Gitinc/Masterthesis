//
// Created by peter on 16.09.21.
//

#ifndef UNTITLED_MACKEYGLASS_H
#define UNTITLED_MACKEYGLASS_H


class mg_parameters{
public:
    double gamma;
    double tau;
    double beta;
    int n;


    mg_parameters(double MGbeta, double MGgamma, int MGn, double MGtau);
};

double deriv(double x, double  x_dot_tau,  mg_parameters MG);

double mackeyglass_rk4(double x_t, double x_t_minus_tau, double deltat, mg_parameters MG);

void mackey(double *X, double *T, int sample_n, mg_parameters MG);

void write_mackey_to_file(int sample_all, int inverval);

    //double mid_step();

    //void write_mackey_to_file(int sample_all);


#endif //UNTITLED_MACKEYGLASS_H
