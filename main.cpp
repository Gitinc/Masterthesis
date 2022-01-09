#include<iostream>
#include<vector>
#include "armadillo"
#include "helper.h"
#include "Mask.h"
#include "MackeyGlass.h"
#include <iomanip>
#include <string>
#include "LorentzAttractor.h"
#include <Narma10.h>
#include <iterator>
#include <random>
#include <chrono>
#include <omp.h>

using namespace std;

class Solution {
public:
    int M;
    int N;
    int tau_max;
    int Nv;
    double c_ring;
    double I_sat;
    double a;
    double K;
    double g;
    double J_0;
    int K_totallength;
    int seed;
    double D;
    int delay;
    string system;
    vector<double> K_values;
    vector<double> x_n;
    vector<vector<double>> K_in;
    vector<vector<double>> K_nin;
    vector<vector<double>> K_out;
    vector<vector<double>> K_training_inputs;


    Solution(int SM, int SN, int Stau, int SNv, double Sc_ring, double SI_sat, double Sa, double SK, double Sg, double SJ_0, int SK_totallength, int Sseed, string Ssystem, double SD, int Sdelay){
        M = SM;
        N = SN;
        tau_max = Stau;
        Nv = SNv;
        c_ring = Sc_ring;
        I_sat = SI_sat;
        a = Sa;
        K = SK;
        g = Sg;
        J_0 = SJ_0;
        D = SD;
        K_totallength = SK_totallength;
        seed = Sseed;
        system = "/home/noah.jaitner/Local/Masterarbeit/Timeseries/" + Ssystem;
        delay = Sdelay;
    }

    void createBasic_Matrix(){
        createK_values();
        createK_in();
        createX_n();
        createK_nin();
        createK_out();
        K_training_inputs = generate_MaskedInputs();
    }

    void createK_in(){
        K_in = createI(N,M);
    }

    void createK_nin(){
        vector<vector<double>> Knin(M, vector<double> (N));

        for (int i_m = 0; i_m < M; i_m++){
            for (int i_n = 0; i_n < N; i_n++){
                if(K_in[i_m][i_n] == 0){
                    Knin[i_m][i_n] = 1;
                }else{
                    Knin[i_m][i_n] = 0;
                }
            }
        }

        K_nin = Knin;
    }

    //K_out works now like in python code, modulos is diffrent in python and c++
    void createK_out(){
        /*std::vector<std::vector<double>> Kout(Nv*tau_max, std::vector<double> (Nv*tau_max));
        std::vector<std::vector<double>> Kout_correct_size(M, std::vector<double> (N));

        int n = 0;
        int m = 0;
        for(int ii = 0; ii < tau_max * Nv; ++ii){
            for(int i = 0; i < tau_max; ++i){
                Kout[ii][(i * Nv) + n] = K_values[((i - m) % tau_max+ tau_max) % tau_max];
                if((((i - m) % tau_max)+tau_max) % tau_max == (tau_max - 1)) {
                    if(((i * Nv) + n) == 0) {
                        Kout[ii][(tau_max * Nv) - 1] = c_ring;
                    }else{
                        Kout[ii][(i * Nv) + n - 1] = c_ring;
                    }
                }
            }
            n +=1;
            if(n == Nv) {
                m += 1;
                n = 0;
            }
        }*/
        K_out = createI(N,M);

        /*for(int i = 0;i<M;i++){
            for(int j = 0;j<N;j++){
                Kout_correct_size[i][j] = Kout[i][j];
            }
        }
        K_out = Kout_correct_size;*/
    }

    void createX_n(){
        vector<double> xn(N, 0);

        for(int i = 0;i<N;i++) xn[i] = 0;

        x_n = xn;
    }

    void createK_values(){
        vector<double> Kvalues(tau_max, 0);

        for(int i = 0;i<Kvalues.size();i++) Kvalues[i] = 0;
        Kvalues[tau_max-1] = 1;
        K_values = Kvalues;
    }

    double nonlinear_G(double x) {
        x = x*exp(-a/(1+x/I_sat));
        //x = a*x/(1+x/I_sat);

        return x;
    }

    vector<double> calculateBackTerm(int m){

        vector<double> Back_Solution(N,0);
        for (int n = 0; n < N; ++n) {
            Back_Solution[n] = K_nin[(m%M+M)%M][n] * (1 - K_out[(m%M+M)%M][n]) * x_n[n];
        }

        return Back_Solution;
    }

    vector<double> calculateFrontTerm(int m, int current_K_training){

        vector<double> Front_Solution(N,0);
        for (int n = 0; n < N; ++n) {
            Front_Solution[n] = K_in[((m-current_K_training)%M+M)%M][n] * nonlinear_G(K * calculateX_out(m-current_K_training) + calculate_currentInput(m,current_K_training));
        }

        return Front_Solution;
    }

    double calculateX_out(int m){
        double x_out_current = 0;
        for(int i = 0;i<N;i++){
            x_out_current += K_out[(m%M+M)%M][i]*x_n[i];
        }

        return x_out_current;
    }

    double calculate_currentInput(int m, int current_K_training){

        double J_k = g*K_training_inputs[current_K_training][m]+J_0;
        return J_k;
    }

    //Has been check, the MaskedInput Matrix is correct
    vector<vector<double>> generate_MaskedInputs(){

        vector<double> X;
        ifstream file(system);
        if (file.is_open()) {
            string line;
            while (getline(file, line)) {
                X.push_back(stod(line));
            }
            file.close();
        }

        Mask Mask_used(Nv,seed);
        Mask Mask_used_delay(Nv,seed+1);

        vector<vector<double>> MaskedInput(K_totallength, vector<double> (Nv));
        int delay_factor = (delay == 0) ? 0 : 1;
        for (int i = 0; i < K_totallength; ++i) {
            for (int j = 0; j < Nv; ++j) {
                MaskedInput[i][j] = X[i]*Mask_used.Mask_vec[j] + delay_factor*X[i+delay]*Mask_used_delay.Mask_vec[j];
            }
        }

        return MaskedInput;
    }

    vector<double> calc_output_forSpecificTrainingStep(int current_K_training){
        vector<double> x_out_Solution(Nv);

        //m needs to be passed in this weird way, due to m in K_out hanging into Nv
        //For example if Nv = 5 and K_out=3x3 after 1 run through nv = 0 again but Index in K_out should be 2
        // Nv: 0 1 2 3 4 0 1 2 3 4...
        // m:  0 1 2 0 1 2 0 1 2 0...
        for (int m = 0; m < Nv; ++m) {
            x_out_Solution[m] = calculateX_out(m-current_K_training);

            vector<double> Back_Solution = calculateBackTerm(m-current_K_training);
            vector<double> Front_Solution = calculateFrontTerm(m,current_K_training);
            for (int n = 0; n < N; ++n) {
                x_n[n] = Front_Solution[n] + Back_Solution[n];

                double mean = 0.0;
                double stddev = x_n[n];
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::default_random_engine generator (seed);
                std::normal_distribution<double> distribution (mean,stddev);

                x_n[n] = x_n[n] + D*distribution(generator);

            }
        }

        return x_out_Solution;
    }

    vector<vector<double>> calcOutputMatrix(){

        vector<vector<double>> OutputMatrix(K_totallength, vector<double> (Nv));

        #pragma omp parallel for
        for (int i = 0; i < K_totallength; ++i) {
                OutputMatrix[i] = calc_output_forSpecificTrainingStep(i);
                OutputMatrix[i].push_back(1);
        }

        /*ofstream ergebnis("Statematrix.txt");
            if (ergebnis.is_open())
            {
                for(int i = 0;i<K_totallength;i++){
                    for(int j = 0; j< Nv;j++){
                        ergebnis << OutputMatrix[i][j] << " ";
                    }
                    ergebnis << "\n";
                }
            }
*/
        return OutputMatrix;
    }
};

/*vector<double> calcMC(Solution used_Model){
        vector<double> random_Uni(75000,0);
        typedef boost::variate_generator<boost::mt19937 &,boost::uniform_real<>> randomGeneratorType;
        boost::mt19937 generator(42u);
                if (input_seed == 0){
                    generator.seed(static_cast<unsigned int>(time(nullptr))); //Seeding the generator
                }
                else{
                    generator.seed(input_seed);
                }
                boost::uniform_real<> uniform(-1.0,1.0);

                randomGeneratorType uni(generator,uniform);

                std::vector<double> random_mask_vector(Nv);

                for (int i = 0; i < Nv; i++)
                {
                    random_mask_vector[i] = uni();
                }
}
*/
vector<double> calculate_System(int K_trainingsteps, int K_teststeps, int ahead_prediciton, int K_buffer, Solution used_Model, bool writeToFile, vector<double>& NRMSE_mean_training, string system_use, int delay){

    system_use = "/home/noah.jaitner/Local/Masterarbeit/Timeseries/" + system_use;
    vector<double> X;
    ifstream file(system_use);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            X.push_back(stod(line));
        }
        file.close();
    }

    vector<double> NRMSE_mean_test(10);
    vector<double> X_Trainingsteps = vector<double>(X.begin() + K_buffer + ahead_prediciton + delay, X.begin() + K_buffer+K_trainingsteps +ahead_prediciton +delay);
    vector<double> X_Teststeps = vector<double>(X.begin() + 2*K_buffer+K_trainingsteps+ahead_prediciton + delay, X.begin()+2*K_buffer+K_teststeps+K_trainingsteps+ahead_prediciton+delay);

    for(int n = 0;n<NRMSE_mean_test.size();n++){

        used_Model.createBasic_Matrix();
        vector<vector<double>> Model_fullStateMatrix = used_Model.calcOutputMatrix();
        arma::mat S_training(K_trainingsteps, used_Model.Nv+1);
        arma::mat S_testing(K_teststeps, used_Model.Nv+1);
        for (int i = 0; i < K_trainingsteps; ++i) {
            for (int j = 0; j < used_Model.Nv+1; ++j) {
                S_training(i,j) = Model_fullStateMatrix[K_buffer+i][j];
            }
        }

        for (int i = 0; i < K_teststeps; ++i) {
            for (int j = 0; j < used_Model.Nv+1; ++j) {
                S_testing(i,j) = Model_fullStateMatrix[2*K_buffer+K_trainingsteps+i][j];
            }
        }

        double lambda = 0.000005;
        arma::dvec w_out = arma::pinv(S_training.t()*S_training+lambda*arma::eye(used_Model.Nv+1,used_Model.Nv+1),0.00000000001)*S_training.t()*arma::conv_to<arma::mat>::from(X_Trainingsteps);
        //arma::dvec w_out = arma::solve(S_training,arma::conv_to<arma::mat>::from(X_Trainingsteps),arma::solve_opts::no_approx);

        arma::mat training_outcome_weights = S_training*w_out;
        arma::mat testing_outcome_weights = S_testing*w_out;
        if(writeToFile){
            ofstream ergebnis("Solutions"+ system_use +".txt");
            if (ergebnis.is_open())
            {
                for(int i = 0; i < K_teststeps; ++i){
                    ergebnis << X_Teststeps[i] << ";"<< testing_outcome_weights(i)<< "\n";
                }
                ergebnis.close();
            }
        }
        NRMSE_mean_test[n] = calc_NRMSE(X_Teststeps, arma::conv_to<vector<double>>::from(testing_outcome_weights));
        NRMSE_mean_training[n] = calc_NRMSE(X_Trainingsteps, arma::conv_to<vector<double>>::from(training_outcome_weights));
    }
    return NRMSE_mean_test;
}

int main (int argc, char *argv[]) {
    //write_mackey_to_file(5000000,100);
    //write_Lorentz_to_file(100000,2);
    //write_Narma10_to_file(50000);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vector<double> NRMSE_mean_training(100);
    vector<double> NRMSE_mean_training_values(10);
    vector<double> NRMSE_mean_test(100);
    vector<double> STD_training(100);
    vector<double> STD_test(100);
    vector<double> current_NRMSE;

    double NRMSE_test;
    double NRMSE_training;

    for (int j = 0; j < 100; ++j) {
        //     Solution(int SM, int SN, int Stau, int SNv, double Sc_ring, double SI_sat, double Sa, double SK, double Sg, double SJ_0, int SK_totallength, int Sseed, string Ssystem, double SD, int delay){
        Solution Test(stoi(argv[1]), stoi(argv[2]), 16, stoi(argv[3]), 1, 0.21910664, 3.00715971, stod(argv[7]), stod(argv[4]), stod(argv[5]), 40000, 0, string(argv[8]), 0.002*j, stoi(argv[9]));
        current_NRMSE = calculate_System(10000,5000,stoi(argv[6]),10000,Test, false, NRMSE_mean_training_values, string(argv[8]),stoi(argv[9]));
        //Solution Test(11,11,16,20,1,1,40,0.01*(j+1),1,0,35000, 0);
        //current_NRMSE = calculate_System(10000,5000,1,10000,Test, true, NRMSE_mean_training_values);
        NRMSE_mean_test[j] = c_mean(current_NRMSE);
        NRMSE_mean_training[j] = c_mean(NRMSE_mean_training_values);
        STD_test[j] = sqrt(Variance(current_NRMSE));
        STD_training[j] = sqrt(Variance(NRMSE_mean_training_values));
    }


    string PathForSolution = string(argv[8]) + string("_MeanWithSTD_M") + string(argv[1]) + string("_N") + string(argv[2]) + string("_Nv") + string(argv[3]) + string("_g") + string(argv[4]) + string("_J0") + string(argv[5]) + string("_Pre") + string(argv[6]) + string("_Delay") + string(argv[9]) + string("_K") + string(argv[7]) + string(".txt");
    ofstream ergebnis(PathForSolution);
    if (ergebnis.is_open())
    {
        for(int i = 0; i < 100; ++i){
            ergebnis << NRMSE_mean_test[i] << " " << STD_test[i] << " " <<  NRMSE_mean_training[i] << " " << STD_training[i] <<"\n";
        }
        ergebnis.close();
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::minutes> (end - begin).count() << "[m]" << std::endl;

    //./Masterarbeit M N Nv g J0 pre D task
    // ko
    // lorenz g<=0.01 sonst keine pinv
}
