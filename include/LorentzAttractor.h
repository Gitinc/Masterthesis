#ifndef LORENTZATTRACTOR_H
#define LORENTZATTRACTOR_H
#include <vector>

void write_Lorentz_to_file(int sample_all, int interval);
void UpdateSolution(std::vector<double> &a);
std::vector<double> deriv(std::vector<double> x);

#endif // LORENTZATTRACTOR_H
