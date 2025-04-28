#include "waveform_Analyser.hh"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TROOT.h> 

#include <regex>

#include "utils.hh"


waveform_Analyser::waveform_Analyser():my_sampling_time(16e-9)
{
    clear_all();
}

waveform_Analyser::waveform_Analyser(std::string file):waveform_Analyser()
{
    this->my_file = file;
    update();
}

waveform_Analyser::waveform_Analyser(const std::vector<double> &adcs):waveform_Analyser()
{
    this->my_adcs = adcs;
}

waveform_Analyser::~waveform_Analyser()
{
}

void waveform_Analyser::set_file(const std::string file)
{
    this->my_file = file;  
}

std::string waveform_Analyser::get_file()
{
    return this->my_file;
}

void waveform_Analyser::set_file_template(const std::string file)
{
    this->my_file_template = file;  
}

std::string waveform_Analyser::get_file_template()
{
    return this->my_file_template;
}

void waveform_Analyser::read_file_template()
{
    this->set_template(this->my_file_template);
}

void waveform_Analyser::set_adcs(const std::vector<double>& adcs)
{
    this->my_adcs = adcs;
}

std::vector<double> waveform_Analyser::get_adcs()
{
    return this->my_adcs;
}

void waveform_Analyser::set_template(const std::vector<double> &template_wf)
{
    this->my_template = template_wf;
}

void waveform_Analyser::set_template(const std::string template_file)
{
    this->my_template.clear();
    std::ifstream file(template_file);  
    
    if (!file.is_open())
    {
        std::cerr << "Error opening the file: " << this->my_file << std::endl;
        return;
    }

    std::string line;
    double value;

    while (std::getline(file, line)) 
    {
        value = std::stod(line);
        this->my_template.push_back(value);
        
    }
    file.close();
}

std::vector<double> waveform_Analyser::get_template()
{
    return this->my_template;
}

void waveform_Analyser::update()
{
    this->my_adcs.clear();
    std::ifstream file(this->my_file);  
    
    if (!file.is_open())
    {
        std::cerr << "Error opening the file: " << this->my_file << std::endl;
        return;
    }

    std::string line;
    double value;

    while (std::getline(file, line)) 
    {
        value = std::stod(line);
        this->my_adcs.push_back(value);
        
    }

    file.close();
}

void waveform_Analyser::conv(const std::vector<double> &template_wf , const int out_size )
{
   this->my_convolution = this->convolve_operation(this->my_adcs,template_wf,out_size);
}

std::vector<double> waveform_Analyser::get_conv()
{
    return this->my_convolution;
}

void waveform_Analyser::ext_conv(const std::vector<double> &wf1, const std::vector<double> &wf2 , const int out_size, const int roll_factor )
{
    this->my_ext_convolution = this->convolve_operation(wf1,wf2,out_size,roll_factor);
}

std::vector<double> waveform_Analyser::get_ext_conv()
{
    return  this->my_ext_convolution;
}

void waveform_Analyser::clear_all()
{
    this->my_file.clear();
    this->my_adcs.clear();
    this->my_convolution.clear();
    this->my_ext_convolution.clear();
    this->my_template.clear();
}

void waveform_Analyser::set_ch(const int ch)
{
    this->my_ch = ch;
}

void waveform_Analyser::set_ch_and_voltage_file()
{
    std::regex pattern("wf_ch(\\d{1,2})_(\\d{1,3})kV\\.txt");
    std::smatch match;

    if (std::regex_search(this->my_file, match, pattern)) {
        int channel = std::stoi(match[1]);
        int voltage = std::stoi(match[2]);

        this->my_ch = channel;
        this->my_voltage = voltage;    
    } 
    else
    {
        std::cout << "File format not recognized!" << std::endl;
    }
}

void waveform_Analyser::set_voltage(const int voltage)
{
    this->my_voltage = voltage;
}

const int waveform_Analyser::get_ch()
{
    return this->my_ch;
}

const int waveform_Analyser::get_voltage()
{
    return this->my_voltage;
}

void waveform_Analyser::set_sampling_time(const double time)
{
    this->my_sampling_time = time;
}

double waveform_Analyser::get_sampling_time()
{
    return this->my_sampling_time;
}

void waveform_Analyser::fit_0(double *params, int n_params)
{
    int size = my_adcs.size();
    double* y_array = &this->my_adcs[0];
    
    std::vector<double> time(size);
    for(int i = 0 ; i < size ; i ++)
    {
        time[i] = i;
    }
    double* x_array = &time[0];

    TGraph *graph = new TGraph(size, x_array, y_array);

    this->my_fit_params_0.clear();
    this->my_fit_params_0 = std::vector<double>(n_params , 0.0);

    for(int i = 0; i< n_params ; i++)
    {
        this->my_fit_params_0[i] = params[i] + 1; //forçar ser diferente na primeira vez, para calcular a 1 convolucao
    }

    TF1 *fitFunc = new TF1("fitFunc_conv", [this](double* x, double* par) { return this->fit_0_function(x, par); }, 0, size, n_params );
    
    for (int i = 0; i < n_params; ++i) 
    {
        fitFunc->SetParameter(i, params[i]);
    }
    this->cont_conv = 0;

    fitFunc->SetParLimits(0,-10000,-0.01);
    fitFunc->SetParLimits(1,-10000,-0.01);
    fitFunc->SetParLimits(2,0.5e-9,100e-9);
    fitFunc->SetParLimits(3,800e-9,2000e-9);
    fitFunc->SetParLimits(5,-80,-40);
    fitFunc->FixParameter(5,params[5]);

    graph->Fit(fitFunc, "0R+"); 

    for (int i = 0; i < n_params; ++i) 
    {
        this->my_fit_params_0[i] = fitFunc->GetParameter(i);
    }

}

const double waveform_Analyser::get_fit_param_0(const int i)
{
    auto n = this->my_fit_params_0.size();
    if (i >= n || i < 0)
    {
        std::cout <<  "Error in the index" << std::endl;
        return -1;
    }
    else
    {
        return this->my_fit_params_0[i];
    }
    return 0.0;
}

const std::vector<double> waveform_Analyser::get_fit_params_0()
{
    return this->my_fit_params_0;
}

void waveform_Analyser::save_fig_fit()
{
    gROOT->SetBatch(kTRUE);
    auto waveform = this->get_adcs();
    auto waveform2 = this->get_ext_conv(); 
    
    int n_points1 = waveform.size();
    int n_points2 = waveform2.size();
    std::vector<double> x1(n_points1);
    std::vector<double> x2(n_points2);
    for (int i = 0; i < n_points1; ++i) 
    {
        x1[i] = i; 
    }
    for (int i = 0; i < n_points2; ++i) 
    {
        x2[i] = i; 
    }

    TCanvas* canvas = new TCanvas("canvas", "Waveform", 800, 600);

    TGraph* graph = new TGraph(n_points2, x2.data(), waveform2.data());
    graph->SetTitle("Waveform;Sample;ADC Value");
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);

    graph->Draw("AL"); 

    TGraph* graph2 = new TGraph(n_points1, x1.data(), waveform.data());
    graph2->SetTitle("Waveform;Sample;ADC Value");
    graph2->SetLineColor(kRed);
    graph2->SetLineWidth(2);

    graph2->Draw("SAME"); 

    canvas->Update();
    canvas->SaveAs(("./figures/waveform_fit_ch" + std::to_string(this->my_ch) + "_" + std::to_string(this->my_voltage) + "kV.png").c_str());

    delete graph;
    delete graph2;
    delete canvas;
    gROOT->SetBatch(kFALSE);

}

std::vector<double> waveform_Analyser::convolve_operation(const std::vector<double> &wf1, const std::vector<double> &wf2, const int out_size, const int roll_factor)
{
    int n1 = wf1.size();
    int n2 = wf2.size();
    
    int n_result = n1 + n2 - 1;
    
    if (out_size != -1)
    {
        n_result = out_size;
    }

    std::vector<double> result(n_result, 0.0);
    
    for (int i = 0; i < n_result; ++i)
    {
        for (int j = 0; j < n1; ++j) 
        {
            if (i - j >= 0 && i - j < n2) 
            {
                result[i] += wf1[j] * wf2[i - j];
            }
        }
    }
    roll_vector(result,roll_factor);
    return result;
}

std::vector<double> waveform_Analyser::LAr_Profile(int n_points, double sample_time, double A1, double A2, double tau1, double tau2, double c)
{
    std::vector<double> s;
    double value; 

    for(int i=0; i < n_points; i++)
    {
        value = A1 * TMath::Exp(-i*sample_time/tau1) + A2 * TMath::Exp(-i*sample_time/tau2) + c ; 
        s.push_back(value);
    }
    return s;
}

double waveform_Analyser::fit_0_function(double *x, double *par)
{
    int n = this->my_template.size();
    int param_size = this->my_fit_params_0.size();

    bool recompute = false;
    for (int i = 0; i < param_size; i++) 
    {
        if (par[i] != this->my_fit_params_0[i]) 
        {
            this->my_fit_params_0[i] = par[i];
            recompute = true;
        }
    }

    if (recompute)
    {
        this->cont_conv++;
        double A1 = par[0];
        double A2 = par[1];
        double tau1 = par[2];
        double tau2 = par[3];
        double c = par[4];
        int roll_factor = (int) par[5];
        auto s = this->LAr_Profile(n, this->my_sampling_time, A1, A2, tau1, tau2 , c);
        this->ext_conv( s , this->my_template , -1 , roll_factor );
        this->my_ext_convolution = std::vector(this->my_ext_convolution.begin(), this->my_ext_convolution.begin() + 1024);
    } 

    int i = static_cast<int>(x[0]);
    if (i >= 0 && i <this->my_ext_convolution.size()) 
    {
        return this->my_ext_convolution[i];
    } 
    else 
    {
        return 0.0;
    }
}
