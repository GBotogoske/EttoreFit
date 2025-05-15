#ifndef UTILS_HH
#define UTILS_HH

#include <vector>

void roll_vector(std::vector<double>& v, int shift); 

std::vector<double> shift_waveform_continuous(const std::vector<double>& input, double shift);

std::vector<double> gaus(double A, double mu, double sigma, int N = 1024);

std::vector<double> norm_gaus(double mu, double sigma, int N = 1024);

#endif
