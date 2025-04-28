#include "waveforms_Analyser.hh"

#include <algorithm>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <TROOT.h> 


namespace fs = std::filesystem;


waveforms_Analyser::waveforms_Analyser() : my_prefix("wf_ch"), my_suffix(".txt") 
{
    this->clear();
}

waveforms_Analyser::~waveforms_Analyser()
{
}

std::string waveforms_Analyser::get_folder()
{
    return this->my_folder_path;
}

void waveforms_Analyser::set_folder(const std::string folder_path)
{
    this->my_folder_path = folder_path;
}

void waveforms_Analyser::update()
{
    std::vector<std::string> matching_files;

    for (const auto& entry : fs::directory_iterator(this->my_folder_path))
    {
        if (entry.is_regular_file()) 
        {
            std::string filename = entry.path().filename().string();
            if (this->starts_with(filename,this->my_prefix) && this->ends_with(filename,this->my_suffix)) 
            {
                matching_files.push_back(entry.path().string());  
            }
        }
    }

    this->wf_Analyser_vector.clear();
    for (int i = 0 ; i < matching_files.size(); i++)
    {
        this->wf_Analyser_vector.push_back(new waveform_Analyser(matching_files[i]));
        this->wf_Analyser_vector[i]->set_ch_and_voltage_file();
        this->wf_Analyser_vector[i]->set_template("template/ch_25_endpoint_112_avg.txt",true); //mudar aqui depois
    }
    std::sort(this->wf_Analyser_vector.begin(), this->wf_Analyser_vector.end(), [this](waveform_Analyser* a, waveform_Analyser* b) {return this->compare_wf(*a, *b);});
    
}

void waveforms_Analyser::clear()
{
    this->my_folder_path = "";
    this->wf_Analyser_vector.clear();
}

std::string waveforms_Analyser::get_suffix()
{
    return this->my_suffix;
}

std::string waveforms_Analyser::get_prefix()
{
    return this->my_prefix;
}

void waveforms_Analyser::set_suffix(const std::string suffix)
{
    this->my_suffix = suffix;
}

void waveforms_Analyser::set_prefix(const std::string prefix)
{
    this->my_prefix = prefix;
}

std::vector<waveform_Analyser *> waveforms_Analyser::get_wf_Analyser_vector()
{
    return this->wf_Analyser_vector;
}

std::vector<waveform_Analyser *> waveforms_Analyser::get_wf_Analyser_vector_by_ch(const int ch)
{
    std::vector<waveform_Analyser *> matching_analyser;
    int n = this->wf_Analyser_vector.size();

    for(int i=0;i<n;i++)
    {
        if(this->wf_Analyser_vector[i]->get_ch() == ch )
        {
            matching_analyser.push_back(this->wf_Analyser_vector[i]);
        }
    }

    return matching_analyser;
}

std::vector<waveform_Analyser *> waveforms_Analyser::get_wf_Analyser_vector_by_voltage(const int voltage)
{
    std::vector<waveform_Analyser *> matching_analyser;
    int n = this->wf_Analyser_vector.size();

    for(int i=0;i<n;i++)
    {
        if(this->wf_Analyser_vector[i]->get_voltage() == voltage )
        {
            matching_analyser.push_back(this->wf_Analyser_vector[i]);
        }
    }

    return matching_analyser;
}

void waveforms_Analyser::plot_slow_comp(std::vector<waveform_Analyser*> analyser)
{
    gROOT->SetBatch(kTRUE);

    int n = analyser.size();
    
    std::vector<double> e_field(n);
    std::vector<double> y(n);

    auto ch=analyser[0]->get_ch();

    for(int i=0; i<n ; i++)
    {
        e_field[i] = analyser[i]->get_voltage();
        y[i] = analyser[i]->get_fit_param_0(3)/1e-9;
    }

    TCanvas* canvas = new TCanvas("canvas", "Slow comp", 800, 600);
    TGraph* graph = new TGraph(n, e_field.data(), y.data());
    graph->SetTitle("Slow comp;Efield [kv_cm];Slow comp [ns]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    graph->GetYaxis()->SetRangeUser(0.98*(*std::min_element(y.begin(), y.end())), 1.02*(*std::max_element(y.begin(), y.end())));

    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(("./figures/" + std::to_string(ch) + "/slow_comp_ch" + std::to_string(ch) + ".png").c_str());

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_fast_comp(std::vector<waveform_Analyser *> analyser)
{
    gROOT->SetBatch(kTRUE);

    int n = analyser.size();
    
    std::vector<double> e_field(n);
    std::vector<double> y(n);

    auto ch=analyser[0]->get_ch();

    for(int i=0; i<n ; i++)
    {
        e_field[i] = analyser[i]->get_voltage();
        y[i] = analyser[i]->get_fit_param_0(2)/1e-9;
    }

    TCanvas* canvas = new TCanvas("canvas", "Fast comp", 800, 600);
    TGraph* graph = new TGraph(n, e_field.data(), y.data());
    graph->SetTitle("Fast comp;Efield [kv_cm];Fast comp [ns]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));   
    graph->GetYaxis()->SetRangeUser(0.9*(*std::min_element(y.begin(), y.end())), 1.1*(*std::max_element(y.begin(), y.end())));
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(("./figures/" + std::to_string(ch) + "/fast_comp_ch" + std::to_string(ch) + ".png").c_str());

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_slow_intensity(std::vector<waveform_Analyser *> analyser)
{
    gROOT->SetBatch(kTRUE);

    int n = analyser.size();
    
    std::vector<double> e_field(n);
    std::vector<double> y(n);

    auto ch=analyser[0]->get_ch();

    for(int i=0; i<n ; i++)
    {
        e_field[i] = analyser[i]->get_voltage();
        y[i] = -analyser[i]->get_fit_param_0(1);
    }

    TCanvas* canvas = new TCanvas("canvas", "Slow I", 800, 600);
    TGraph* graph = new TGraph(n, e_field.data(), y.data());
    graph->SetTitle("Slow intensity ; Efield [kv_cm] ; Slow intensity");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    graph->GetYaxis()->SetRangeUser(0.9*(*std::min_element(y.begin(), y.end())), 1.1*(*std::max_element(y.begin(), y.end())));
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(("./figures/" + std::to_string(ch) + "/slow_intensity_ch" + std::to_string(ch) + ".png").c_str());

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_fast_intensity(std::vector<waveform_Analyser *> analyser)
{
    gROOT->SetBatch(kTRUE);

    int n = analyser.size();
    
    std::vector<double> e_field(n);
    std::vector<double> y(n);

    auto ch=analyser[0]->get_ch();

    for(int i=0; i<n ; i++)
    {
        e_field[i] = analyser[i]->get_voltage();
        y[i] = -analyser[i]->get_fit_param_0(0);
    }

    TCanvas* canvas = new TCanvas("canvas", "Fast I", 800, 600);
    TGraph* graph = new TGraph(n, e_field.data(), y.data());
    graph->SetTitle("Fast intensity ; Efield [kv_cm] ; Fast intensity ");
   
    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end())); 
    graph->GetYaxis()->SetRangeUser(0.9*(*std::min_element(y.begin(), y.end())), 1.1*(*std::max_element(y.begin(), y.end())));
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(("./figures/" + std::to_string(ch) + "/fast_intensity_ch" + std::to_string(ch) + ".png").c_str());

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_all(const int ch)
{

    std::string folder_path = "figures/" + std::to_string(ch);

    if (!fs::exists(folder_path)) 
    {
        fs::create_directory(folder_path);
    }

    std::vector<waveform_Analyser *> matching_analyser = this->get_wf_Analyser_vector_by_ch(ch);
    this->plot_slow_comp(matching_analyser);
    this->plot_fast_comp(matching_analyser);
    this->plot_slow_intensity(matching_analyser);
    this->plot_fast_intensity(matching_analyser);
    
}

void waveforms_Analyser::plot_all_channels(const int n_ch)
{
    for(int i=0; i < n_ch ; i++)
    {
        this->plot_all(i);
    }
}

void waveforms_Analyser::fit_channel(const int ch, double *par, int n_par)
{
    std::vector<waveform_Analyser *> matching_analyser = this->get_wf_Analyser_vector_by_ch(ch);
    int n = matching_analyser.size();

    for(int i=0; i<n ; i++)
    {
        matching_analyser[i]->fit_0(par,n_par);
        matching_analyser[i]->save_fig_fit();
    }
}

void waveforms_Analyser::fit_all_channels(const int n_ch, double *par, int n_par)
{
    for(int i=0; i < n_ch ; i++)
    {
        this->fit_channel(i,par,n_par);
    }
}

bool waveforms_Analyser::starts_with(const std::string &str, const std::string &prefix)
{
    return str.substr(0, prefix.length()) == prefix;
}

bool waveforms_Analyser::ends_with(const std::string &str, const std::string &suffix)
{
    if (str.length() < suffix.length()) 
    {
        return false;
    }
    return str.substr(str.length() - suffix.length()) == suffix;
}

bool waveforms_Analyser::compare_wf(waveform_Analyser& wf1, waveform_Analyser& wf2)
{
    auto a = wf1.get_ch();
    auto b = wf2.get_ch();
    auto a_v = wf1.get_voltage();
    auto b_v = wf2.get_voltage();
    if (a != b)
        return a < b; // ordena por canal
    else
        return a_v < b_v;
}
