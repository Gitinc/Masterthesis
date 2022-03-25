#include "helper.h"
#include<vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

std::vector<double> multiplyMatrixVector(std::vector<std::vector<double>> A,std::vector<double> u){
    std::vector<double> AmultiU(u.size(),0);

    for(int i = 0; i<A[0].size();i++){
        for(int j = 0; j<u.size();j++){
            AmultiU[i] += A[i][j]*u[j];
        }
    }
    return AmultiU;
}

std::vector<double> multiplyScalarVector(double lambda,std::vector<double> u){
    std::vector<double> UmultiLambda(u.size(),0);

        for(int i = 0; i<u.size();i++){
            UmultiLambda[i] += lambda*u[i];
        }

    return UmultiLambda;
}

std::vector<std::vector<double>> multiplyScalarMatrix(double lambda,std::vector<std::vector<double>> A){
    std::vector<double> row(A.size(), 0);
    std::vector<std::vector<double>> Amultilambda(A[0].size(), row);

    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[0].size();j++)
        {
            Amultilambda[i][j]=lambda*A[i][j];
        }
    }

    return Amultilambda;
}

std::vector<std::vector<double>> multiplyMatrix(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B){
    if ( A.size() != B.size() || A[0].size()!=B[0].size() ) {
        throw std::invalid_argument( "Matrixes are not of same Size" );
    }

    std::vector<double> row(A.size(), 0);

    std::vector<std::vector<double>> AmultiB(A[0].size(), row);

    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[0].size();j++)
        {
            AmultiB[i][j]=0;
            for(int k=0;k<A[0].size();k++)
            {
                AmultiB[i][j]+=A[i][k]*B[k][j];
            }
        }
    }

    return AmultiB;
}

std::vector<std::vector<double>> substractMatrix(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B){
    if ( A.size() != B.size() || A[0].size()!=B[0].size() ) {
        throw std::invalid_argument( "Matrixes are not of same Size" );
    }

    std::vector<double> row(A.size(), 0);

    std::vector<std::vector<double>> AminusB(A[0].size(), row);

    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[0].size();j++)
        {
            AminusB[i][j] = A[i][j]-B[i][j];
        }
    }

    return AminusB;
}

std::vector<std::vector<double>> createI(int n, int m){
    std::vector<std::vector<double>> I(m, std::vector<double> (n));

    for (int i = 0; i < m; ++i) {
        for(int j = 0;j<n;++j) {
            if(i==j) I[i][j] = 1;
        }
    }

    return I;
}

std::vector<std::vector<double>> ShiftMatrix(std::vector<std::vector<double>> matrix)
{
    for (auto &row: matrix) // move columns to the left
    {
        rotate(row.begin(), row.begin() + 1, row.end());
    }
    return matrix;
}

double c_mean(std::vector<double> input){
    double sum = 0;
    for (int i = 0; i < input.size(); ++i) {
        sum += input[i];
    }
    return sum/input.size();
}

double Variance(std::vector<double> o){
    double var = 0;
    double mean = c_mean(o);
    for (int i = 0; i < o.size(); ++i) {
        var += pow(mean-o[i],2);
    }
    var /= o.size();
    return var;
}

double calc_NRMSE(std::vector<double> o_expected, std::vector<double> o_calculated){
    double NRMSE = 0;
    for (int i = 0; i < o_expected.size(); ++i) {
        NRMSE += pow(o_expected[i]-o_calculated[i],2);
    }
    NRMSE = sqrt(NRMSE/(o_expected.size()*Variance(o_expected)));

    return NRMSE;
}
