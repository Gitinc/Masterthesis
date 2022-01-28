#include "Narma10.h"
#include<boost/random.hpp>
#include <sys/time.h>
#include "iostream"
#include <fstream>
#include <iomanip>

Narma10::Narma10(double N10alpha, double N10beta, double N10gamma, double N10delta)
{
    //ctor
    alpha = N10alpha;
    beta = N10beta;
    gamma = N10gamma;
    delta = N10delta;
}

void Narma10::initalize_RandomNumbers(){
    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);
    time_t msecs_time = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);
    typedef boost::variate_generator<boost::mt19937 &,boost::uniform_real<>> randomGeneratorType;
    boost::mt19937 generator(42u);
    generator.seed(static_cast<unsigned int>(msecs_time));
    boost::uniform_real<> uniform(0.0,0.5);
    randomGeneratorType uni(generator,uniform);
    for (int i = 0; i < sizeof(u_history)/sizeof(u_history[0]); ++i)
    {
        u_history[i] = uni();
    }
}

void Narma10::updateY_k(int seed_rand){
    struct timeval time_now{};
    gettimeofday(&time_now, nullptr);
    time_t msecs_time = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);
    typedef boost::variate_generator<boost::mt19937 &,boost::uniform_real<>> randomGeneratorType;
    boost::mt19937 generator(42u);
    generator.seed(static_cast<unsigned int>(msecs_time+seed_rand));
    boost::uniform_real<> uniform(0.0,0.5);
    randomGeneratorType uni(generator,uniform);

    double tmp = 0.0;
    double sum_OfHistory = 0.0;

    for(int i = 0; i<sizeof(y_history)/sizeof(y_history[0]);i++){
        sum_OfHistory += y_history[i];
    }

    tmp = alpha*y_history[9] + beta*y_history[9]*sum_OfHistory+gamma*u_history[9]*u_history[0]+delta;

    if(tmp > 1) tmp = 1;

    for(int i = 0; i<sizeof(y_history)/sizeof(y_history[0])-1;i++){
        y_history[i] = y_history[i+1];
        u_history[i] = u_history[i+1];
    }
    y_history[9] = tmp;
    u_history[9] = uni();

}

void write_Narma10_to_file(int sample_all){
    Narma10 TimeSeries(0.3,0.05,1.5,0.1);
    double X[sample_all];
    double U[sample_all];
    for(int i = 0;i< sample_all;i++){
        TimeSeries.updateY_k(i);
        X[i] = TimeSeries.y_history[9];
        U[i] = TimeSeries.u_history[9];
    }

    std::ofstream myfile ("narma.dat");
    if (myfile.is_open())
    {
        for(int i = 0; i < sample_all; ++i){
            myfile << std::fixed << std::setprecision(15) << X[i] << "\n" ;
        }
        myfile.close();
    }
    else  throw std::invalid_argument( "File not open" );
    std::ofstream file ("narma_Input.dat");
    if (file.is_open())
    {
        for(int i = 0; i < sample_all; ++i){
            file << std::fixed << std::setprecision(15) << U[i] << "\n" ;
        }
        file.close();
    }
    else  throw std::invalid_argument( "File not open" );
}
