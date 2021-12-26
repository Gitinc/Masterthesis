//
// Created by peter on 18.09.21.
//

#ifndef UNTITLED_MASK_H
#define UNTITLED_MASK_H

class Mask{
public:
    int Nv;
    int input_seed;
    std::vector<double> Mask_vec;

    Mask(int size, int Mseed);

    void generateMask();
};

#endif //UNTITLED_MASK_H
