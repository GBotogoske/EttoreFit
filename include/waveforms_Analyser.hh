#ifndef WAVEFORMS_ANALYSER_HH
#define WAVEFORMS_ANALYSER_HH

#include "waveform_Analyser.hh"
#include <vector>
#include <string>
#include <filesystem>

class waveforms_Analyser
{

    public:
        waveforms_Analyser();
        ~waveforms_Analyser();
        std::string get_folder();
        void set_folder(const std::string folder_path);
        void set_template_file(const std::string folder_path);
        void update();
        void clear();

        std::string get_suffix();
        std::string get_prefix();
        void set_suffix(const std::string suffix);
        void set_prefix(const std::string prefix);

        std::vector<waveform_Analyser*> get_wf_Analyser_vector();
        std::vector<waveform_Analyser*> get_wf_Analyser_vector_by_ch(const int ch);
        std::vector<waveform_Analyser*> get_wf_Analyser_vector_by_voltage(const int voltage);

        std::vector<double> get_par_fit_ch(const int ch, const int par);
        std::vector<double> get_slow_comp_vector(const int ch);
        std::vector<double> get_fast_comp_vector(const int ch);
        std::vector<double> get_slow_int_vector(const int ch);
        std::vector<double> get_fast_int_vector(const int ch);
        std::vector<double> get_ly_vector(const int ch);


        void plot_slow_comp_avg(std::vector<int> list_ch);
        void plot_fast_comp_avg(std::vector<int> list_ch);
        void plot_slow_int_avg(std::vector<int> list_ch);
        void plot_fast_int_avg(std::vector<int> list_ch);
        void plot_all_avg(std::vector<int> list_ch);
        void plot_slow_comp(std::vector<waveform_Analyser*> analyser);
        void plot_fast_comp(std::vector<waveform_Analyser*> analyser);
        void plot_slow_intensity(std::vector<waveform_Analyser*> analyser);
        void plot_fast_intensity(std::vector<waveform_Analyser*> analyser);
        void plot_all(const int ch);
        void plot_all_channels(const int n_ch);
        void fit_channel(const int ch, double* par, int n_par);
        void fit_all_channels(const int n_ch,double* par, int n_par);

        void plot_ly_ch(const int ch);
        void plot_ly_avg(std::vector<int> list_ch);

    private:
        std::vector<waveform_Analyser*> wf_Analyser_vector;
        std::string my_folder_path;
        std::string my_template_path;

        std::string my_prefix;
        std::string my_suffix;
        bool starts_with(const std::string& str, const std::string& prefix);
        bool ends_with(const std::string& str, const std::string& suffix);
        bool compare_wf(waveform_Analyser& wf1, waveform_Analyser& wf2);
        
};


#endif