#ifndef WAVEFORM_ANALYSER_HH
#define WAVEFORM_ANALYSER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class waveform_Analyser
{
    public:
        waveform_Analyser();
        waveform_Analyser(std::string file);
        waveform_Analyser(const std::vector <double>& adcs);
        ~waveform_Analyser();

        void set_file(const std::string file);
        std::string get_file();

        void set_file_template(const std::string file);
        std::string get_file_template();
        void read_file_template();

        void set_adcs(const std::vector <double>& adcs);
        std::vector <double> get_adcs();

        void set_template(const std::vector <double>& template_wf);
        void set_template(const std::string template_file, bool negative = true);
        std::vector <double> get_template();

        void update();

        void conv(const std::vector <double>& template_wf , const int out_size = -1);
        std::vector <double> get_conv();
        void ext_conv(const std::vector <double>& wf1, const std::vector <double>& wf2 , const int out_size = -1, const int roll_factor = 0);
        std::vector <double> get_ext_conv();

        void clear_all();

        void set_ch(const int ch);
        virtual void set_ch_and_voltage_file();
        void set_voltage(const int voltage);
        const int get_ch();
        const int get_voltage();

        void set_sampling_time(const double time);
        double get_sampling_time();

        void fit_0(double *params, int n_params); // do the fit using the ext conv
        void fit_0_discrete_p5_scan(double *params, int n_params);

        const double get_fit_param_0(const int i);
        const std::vector<double> get_fit_params_0();

        const double calc_light_yield();

        void save_fig_fit();

        int cont_conv = 0;

    private:
        std::string my_file;
        std::string my_file_template;
        std::vector <double> my_adcs;
        std::vector <double> my_template;
        std::vector <double> my_convolution;
        std::vector <double> my_ext_convolution;
        std::vector <double> convolve_operation(const std::vector <double>& wf1, const std::vector <double>& wf2, const int out_size = -1, const int roll_factor = 0);

        int my_ch;
        int my_voltage;

        double my_sampling_time;

        std::vector<double> LAr_Profile(int n_points=1024, double sample_time = 16e-9, double A1=1, double A2=1, double tau1=2e-9, double tau2=1600e-9, double c = 0);
        double fit_0_function(double *x, double *par);
        std::vector<double> my_fit_params_0;

        
};

#endif