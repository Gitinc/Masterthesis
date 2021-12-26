#include "vector"
#include "Mask.h"
#include<boost/random.hpp>
#include <sys/time.h>
#include <iostream>

    Mask::Mask(int size, int Mseed){
        Nv = size;
        input_seed = Mseed;
        generateMask();
    }

    void Mask::generateMask(){
        struct timeval time_now{};
        gettimeofday(&time_now, nullptr);
        time_t msecs_time = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);

        std::vector<double> generated_Mask(Nv,0);
        typedef boost::variate_generator<boost::mt19937 &,boost::uniform_real<>> randomGeneratorType;
        boost::mt19937 generator(42u);
                if (input_seed == 0){
                    generator.seed(static_cast<unsigned int>(msecs_time)); //Seeding the generator
                }
                else{
                    generator.seed(input_seed);
                }
                boost::uniform_real<> uniform(0,1.0);

                randomGeneratorType uni(generator,uniform);

                std::vector<double> random_mask_vector(Nv);

                for (int i = 0; i < Nv; i++)
                {
                    random_mask_vector[i] = uni();
                }
                Mask_vec = random_mask_vector;
            }


