#include "utils.hh"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numbers> 

void roll_vector(std::vector<double>& v, int shift)
{
    int n = v.size();
    if (n == 0) return;

    shift = ((shift % n) + n) % n;

    std::rotate(v.begin(), v.end() - shift, v.end());
}

std::vector<double> gaus(double A, double mu, double sigma, int N)
{
    std::vector<double> x(N);
    std::vector<double> my_gaus(N);
    A=std::abs(A);
    // Preencher x com valores igualmente espaçados em torno de mu
    double step = 1.0;  // passo entre pontos (ajuste como necessário)
    double start = mu - (N / 2) * step;

    for(int i = 0; i < N; i++)
    {
        x[i] = start + i * step;
        my_gaus[i] = A * std::exp(-std::pow(x[i] - mu, 2) / (2 * sigma * sigma));
    }

    return my_gaus;
}

std::vector<double> norm_gaus(double mu, double sigma, int N)
{
    sigma=std::abs(sigma);
    return gaus(1/(sigma*std::sqrt(2*M_PI)),mu,sigma,N);
}

std::vector<double> shift_waveform_continuous(const std::vector<double>& input, double shift)
{
    int n = input.size();
    std::vector<double> output(n, 0.0);

    // Corrige shift negativo ou > n
    while (shift < 0) shift += n;
    while (shift >= n) shift -= n;

    for (int i = 0; i < n; ++i)
    {
        double src_idx = i + shift;
        int idx_low = static_cast<int>(std::floor(src_idx)) % n;
        int idx_high = (idx_low + 1) % n;
        double frac = src_idx - std::floor(src_idx);

        // Interpolação linear circular
        output[i] = (1.0 - frac) * input[idx_low] + frac * input[idx_high];
    }

    return output;
}
