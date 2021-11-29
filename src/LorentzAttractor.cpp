#include "LorentzAttractor.h"
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include "iostream"
#include <iomanip>

std::vector<double> deriv(std::vector<double> x){
    std::vector<double> abl(3,0);
    double rho = 28;
    double beta = 8/3;
    double sigma = 10;

    abl[0] = sigma * (x[1]-x[0]);
    abl[1] = x[0]*(rho-x[2])-x[1];
    abl[2] = (x[1]*x[0])-(beta * x[2]);

    return abl;
}

void UpdateSolution(std::vector<double> &a){

    double h = 0.01;
    std::vector<double> k1(a.size());
	std::vector<double> k2(a.size());
	std::vector<double> k3(a.size());
	std::vector<double> k4(a.size());
	std::vector<double> tmp(a.size());

	for(int i = 0;i<a.size();i++){
        k1[i] = h * deriv(a)[i];

        for(int j = 0;j<a.size();j++) tmp[j] = a[j] + 0.5*k1[j];

        k2[i] = h * deriv(tmp)[i];

        for(int j = 0;j<a.size();j++) tmp[j] = a[j] + 0.5*k2[j];

        k3[i] = h * deriv(tmp)[i];

        for(int j = 0;j<a.size();j++) tmp[j] = a[j] + k3[j];

        k4[i] = h * deriv(tmp)[i];

        a[i] = a[i] + (1.0/6.0) * (k1[i] + (2.0 * k2[i]) + (2.0 * k3[i]) + k4[i]);
	}
}

void write_Lorentz_to_file(int sample_all, int interval){
    std::vector<double> a(3);
    a[0] = 1;
    a[1] = 0.5;
    a[2] = -0.3;
    double X[(int) round(sample_all/interval)];
    for (int i = 0; i < (int) round(sample_all/interval); ++i) X[i] = 0.0;
    double Y[(int) round(sample_all/interval)];
    for (int i = 0; i < (int) round(sample_all/interval); ++i) Y[i] = 0.0;
    double Z[(int) round(sample_all/interval)];
    for (int i = 0; i < (int) round(sample_all/interval); ++i) Z[i] = 0.0;

    for(int i = 0; i < sample_all; i++){
        if(i%interval == 0) X[(int) round(i/interval)] = a[0];
        if(i%interval == 0) Y[(int) round(i/interval)] = a[1];
        if(i%interval == 0) Z[(int) round(i/interval)] = a[2];
		UpdateSolution(a);
    }

    std::ofstream myfile ("lorenz.dat");
    if (myfile.is_open())
    {
        for(int i = 0; i < (int) round(sample_all/interval); ++i){
            myfile << std::fixed << std::setprecision(15) << X[i] << "\n" ;
        }
        myfile.close();
    }
    else  throw std::invalid_argument("File not open");
}
