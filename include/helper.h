//
// Created by peter on 11.09.21.
//

#ifndef UNTITLED_HELPER_H
#define UNTITLED_HELPER_H
#include<vector>
#include <stdexcept>

std::vector<double> multiplyMatrixVector(std::vector<std::vector<double>> A,std::vector<double> u);
std::vector<std::vector<double>> multiplyMatrix(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B);
std::vector<std::vector<double>> substractMatrix(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B);
std::vector<std::vector<double>> createI(int n, int m);
std::vector<double> multiplyScalarVector(double lambda,std::vector<double> u);
std::vector<std::vector<double>> ShiftMatrix(std::vector<std::vector<double>> matrix);
double c_mean(std::vector<double> input);
double Variance(std::vector<double> o);
double calc_NRMSE(std::vector<double> o_expected, std::vector<double> o_calculated);


#endif //UNTITLED_HELPER_H