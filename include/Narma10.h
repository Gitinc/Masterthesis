#ifndef NARMA10_H
#define NARMA10_H


class Narma10
{
    public:
        Narma10(double N10alpha, double N10beta, double N10gamma, double N10delta);
        double alpha, beta, gamma, delta;
        double u_history [10] = {0,0,0,0,0,0,0,0,0,0};
        void updateY_k(int seed_rand);
        void initalize_RandomNumbers();
        double y_history [10] = {0,0,0,0,0,0,0,0,0,0};

    protected:

    private:
};

void write_Narma10_to_file(int sample_all);


#endif // NARMA10_H
