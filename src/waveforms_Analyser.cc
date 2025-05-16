#include "waveforms_Analyser.hh"

#include <algorithm>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <TROOT.h> 
#include "TGraphErrors.h"


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

void waveforms_Analyser::set_template_file(const std::string folder_path)
{
    this->my_template_path = folder_path;
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

    std::ifstream file_template(this->my_template_path);
    std::vector<std::string> template_files;
    std::string lines;

    if (!file_template.is_open()) {
        std::cerr << "Template file not opened" << std::endl;
        return;
    }

    while (std::getline(file_template, lines)) {
        template_files.push_back(lines);
    }

    file_template.close();

    this->wf_Analyser_vector.clear();
    for (int i = 0 ; i < matching_files.size(); i++)
    {
        this->wf_Analyser_vector.push_back(new waveform_Analyser(matching_files[i]));
        this->wf_Analyser_vector[i]->set_ch_and_voltage_file();
        int ch = wf_Analyser_vector[i]->get_ch();
        this->wf_Analyser_vector[i]->set_template(template_files[ch],true); 
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

std::vector<double> waveforms_Analyser::get_par_fit_ch(const int ch, const int par)
{
    std::vector<waveform_Analyser *> matching_analyser = this->get_wf_Analyser_vector_by_ch(ch);
    int n = matching_analyser.size();
    
    std::vector<double> y(n);

    for(int i=0; i<n ; i++)
    {
        y[i] = matching_analyser[i]->get_fit_param_0(par);
    }

    return y;
}

std::vector<double> waveforms_Analyser::get_slow_comp_vector(const int ch)
{
    auto y = this->get_par_fit_ch(ch,3);
    int n = y.size();
    for(int i=0; i<n ; i++)
    {
        y[i] = y[i]/1e-9;
    }
    return y;
}

std::vector<double> waveforms_Analyser::get_fast_comp_vector(const int ch)
{
    auto y = this->get_par_fit_ch(ch,2);
    int n = y.size();
    for(int i=0; i<n ; i++)
    {
        y[i] = y[i]/1e-9;
    }
    return y;
}

std::vector<double> waveforms_Analyser::get_slow_int_vector(const int ch)
{
    auto y = this->get_par_fit_ch(ch,1);
    int n = y.size();
    for(int i=0; i<n ; i++)
    {
        y[i] = -y[i];
    }
    return y;
}

std::vector<double> waveforms_Analyser::get_fast_int_vector(const int ch)
{
    auto y = this->get_par_fit_ch(ch,0);
    int n = y.size();
    for(int i=0; i<n ; i++)
    {
        y[i] = -y[i];
    }
    return y;
}

std::vector<double> waveforms_Analyser::get_ly_vector(const int ch)
{
    std::vector<waveform_Analyser *> matching_analyser = this->get_wf_Analyser_vector_by_ch(ch);

    std::vector <double> y;

    for(waveform_Analyser* my_analyser: matching_analyser)
    {   
        y.push_back(my_analyser->calc_light_yield());
    }

    int n = y.size();

    double norm=y[0];
    for(int i=0;i<n;i++)
    {
        y[i]=y[i]/norm;
    }

    return y;
}

void waveforms_Analyser::plot_slow_comp_avg(std::vector<int> list_ch)
{
    std::vector<std::vector<double>> vector;
    for(const int& this_ch : list_ch)
    {
        vector.push_back(this->get_slow_comp_vector(this_ch));
    }

    std::size_t n = vector.size(); //number of channels
    std::size_t dim = vector[0].size(); //nummber of voltage
    std::vector<double> mean(dim, 0.0);
    std::vector<double> stddev(dim, 0.0);

   
    for (const std::vector<double>& this_vector : vector) 
    {
        for (std::size_t i = 0; i < dim; ++i) 
        {
            mean[i] += this_vector[i];
        }
    }
    for (std::size_t i = 0; i < dim; ++i)
    {
        mean[i] /= n;
    }

    // Calcula o desvio padrão
    for (const std::vector<double>& this_vector : vector) 
    {
        for (std::size_t i = 0; i < dim; ++i)
        {
            double diff = this_vector[i] - mean[i];
            stddev[i] += diff * diff;
        }
    }
    for (std::size_t i = 0; i < dim; ++i)
    {
        stddev[i] = std::sqrt(stddev[i] / n); // para amostral: dividir por (n - 1)
    }
  
   
    auto aux_analyser = this->get_wf_Analyser_vector_by_ch(list_ch[0]);
    std::vector<double> e_field(dim);
    for(int i=0; i<dim ; i++)
    {
        e_field[i] = aux_analyser[i]->get_voltage();  
    }

    std::vector<double> ex(dim, 0.0);     

    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("canvas", "Slow comp Avg", 800, 600);
    TGraphErrors* graph = new TGraphErrors(dim, e_field.data(), mean.data(), ex.data(), stddev.data());
    graph->SetTitle("Slow comp; Efield [kv_cm]; Slow comp [ns]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    double y_min = 1e9, y_max = -1e9;
    for (int i = 0; i < n; ++i)
    {
        double y1 = mean[i] - stddev[i];
        double y2 = mean[i] + stddev[i];
        if (y1 < y_min) y_min = y1;
        if (y2 > y_max) y_max = y2;
    }
   //graph->GetYaxis()->SetRangeUser(y_min , y_max );


    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(std::string("./figures/avg/slow_comp.png").c_str()); 

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_fast_comp_avg(std::vector<int> list_ch)
{
    std::vector<std::vector<double>> vector;
    for(const int& this_ch : list_ch)
    {
        vector.push_back(this->get_fast_comp_vector(this_ch));
    }

    std::size_t n = vector.size();
    std::size_t dim = vector[0].size();
    std::vector<double> mean(dim, 0.0);
    std::vector<double> stddev(dim, 0.0);

   
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            mean[i] += this_vector[i];
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        mean[i] /= n;
    }

    // Calcula o desvio padrão
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            double diff = this_vector[i] - mean[i];
            stddev[i] += diff * diff;
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        stddev[i] = std::sqrt(stddev[i] / n); // para amostral: dividir por (n - 1)
    }
  
    auto aux_analyser = this->get_wf_Analyser_vector_by_ch(list_ch[0]);
    std::vector<double> e_field(dim);
    for(int i=0; i<dim ; i++)
    {
        e_field[i] = aux_analyser[i]->get_voltage();  
    }

    std::vector<double> ex(dim, 0.0);     

    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("canvas", "Fast Comp Average", 800, 600);
    TGraphErrors* graph = new TGraphErrors(dim, e_field.data(), mean.data(), ex.data(), stddev.data());
    graph->SetTitle("Fast comp;Efield [kv_cm];Fast comp [ns]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    double y_min = 1e9, y_max = -1e9;
    for (int i = 0; i < n; ++i) {
        double y1 = mean[i] - stddev[i];
        double y2 = mean[i] + stddev[i];
        if (y1 < y_min) y_min = y1;
        if (y2 > y_max) y_max = y2;
    }
   //graph->GetYaxis()->SetRangeUser(y_min , y_max );

    

    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(std::string("./figures/avg/fast_comp.png").c_str()); 

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);

}

void waveforms_Analyser::plot_slow_int_avg(std::vector<int> list_ch)
{
    std::vector<std::vector<double>> vector;
    for(const int& this_ch : list_ch)
    {
        vector.push_back(this->get_slow_int_vector(this_ch));
    }

    std::size_t n = vector.size();
    std::size_t dim = vector[0].size();
    std::vector<double> mean(dim, 0.0);
    std::vector<double> stddev(dim, 0.0);

   
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            mean[i] += this_vector[i];
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        mean[i] /= n;
    }

    // Calcula o desvio padrão
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            double diff = this_vector[i] - mean[i];
            stddev[i] += diff * diff;
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        stddev[i] = std::sqrt(stddev[i] / n); // para amostral: dividir por (n - 1)
    }
  
    auto aux_analyser = this->get_wf_Analyser_vector_by_ch(list_ch[0]);
    std::vector<double> e_field(dim);
    for(int i=0; i<dim ; i++)
    {
        e_field[i] = aux_analyser[i]->get_voltage();  
    }

    std::vector<double> ex(dim, 0.0);   

    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("canvas", "Slow Intensity Avg", 800, 600);
    TGraphErrors* graph = new TGraphErrors(dim, e_field.data(), mean.data(), ex.data(), stddev.data());
    graph->SetTitle("Slow intensity;Efield [kv_cm];Slow Intensity [A.U.]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    double y_min = 1e9, y_max = -1e9;
    for (int i = 0; i < n; ++i) {
        double y1 = mean[i] - stddev[i];
        double y2 = mean[i] + stddev[i];
        if (y1 < y_min) y_min = y1;
        if (y2 > y_max) y_max = y2;
    }
   //graph->GetYaxis()->SetRangeUser(y_min , y_max );

    
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(std::string("./figures/avg/slow_int.png").c_str()); 

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);

}

void waveforms_Analyser::plot_fast_int_avg(std::vector<int> list_ch)
{
    std::vector<std::vector<double>> vector;
    for(const int& this_ch : list_ch)
    {
        vector.push_back(this->get_fast_int_vector(this_ch));
    }

    std::size_t n = vector.size();
    std::size_t dim = vector[0].size();
    std::vector<double> mean(dim, 0.0);
    std::vector<double> stddev(dim, 0.0);

   
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            mean[i] += this_vector[i];
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        mean[i] /= n;
    }

    // Calcula o desvio padrão
    for (const std::vector<double>& this_vector : vector) {
        for (std::size_t i = 0; i < dim; ++i) {
            double diff = this_vector[i] - mean[i];
            stddev[i] += diff * diff;
        }
    }
    for (std::size_t i = 0; i < dim; ++i) {
        stddev[i] = std::sqrt(stddev[i] / n); // para amostral: dividir por (n - 1)
    }
  
    auto aux_analyser = this->get_wf_Analyser_vector_by_ch(list_ch[0]);
    std::vector<double> e_field(dim);
    for(int i=0; i<dim ; i++)
    {
        e_field[i] = aux_analyser[i]->get_voltage();  
    }

    std::vector<double> ex(dim, 0.0);     

    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("canvas", "Fast intensity Avg", 800, 600);
    TGraphErrors* graph = new TGraphErrors(dim, e_field.data(), mean.data(), ex.data(), stddev.data());
    graph->SetTitle("Fast Int;Efield [kv_cm];Fast Int [A.U.]");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    double y_min = 1e9, y_max = -1e9;
    for (int i = 0; i < n; ++i) {
        double y1 = mean[i] - stddev[i];
        double y2 = mean[i] + stddev[i];
        if (y1 < y_min) y_min = y1;
        if (y2 > y_max) y_max = y2;
    }
   //graph->GetYaxis()->SetRangeUser(y_min , y_max );

    
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(std::string("./figures/avg/fast_int.png").c_str()); 

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);

}

void waveforms_Analyser::plot_all_avg(std::vector<int> list_ch)
{
    this->plot_slow_comp_avg(list_ch);
    this->plot_fast_comp_avg(list_ch);
    this->plot_slow_int_avg(list_ch);
    this->plot_fast_int_avg(list_ch);
    this->plot_ly_avg(list_ch);
}

void waveforms_Analyser::plot_slow_comp(std::vector<waveform_Analyser *> analyser)
{
    gROOT->SetBatch(kTRUE);

    int n = analyser.size();
    
    std::vector<double> e_field(n);
    std::vector<double> y(n);
    std::vector<double> ex(n, 0.0);         // Erros no eixo X (zero)
    std::vector<double> ey(n, 0.0);       // Erros no eixo Y (100 ns)

    auto ch=analyser[0]->get_ch();

    for(int i=0; i<n ; i++)
    {
        e_field[i] = analyser[i]->get_voltage();
        y[i] = analyser[i]->get_fit_param_0(3)/1e-9;
    }

    TCanvas* canvas = new TCanvas("canvas", "Slow comp", 800, 600);
    TGraphErrors* graph = new TGraphErrors(n, e_field.data(), y.data(), ex.data(), ey.data());
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
    this->plot_ly_ch(ch);
    
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
        this->plot_all(i);
    }
}

void waveforms_Analyser::plot_ly_ch(const int ch)
{
    std::vector<waveform_Analyser *> matching_analyser = this->get_wf_Analyser_vector_by_ch(ch);

    std::vector <double> y;
    std::vector <double> e_field;

    int j=0;
    for(waveform_Analyser* my_analyser: matching_analyser)
    {   
        y.push_back(my_analyser->calc_light_yield());
        e_field.push_back(my_analyser->get_voltage());
        std::cout <<  y[j] << " --- " << e_field[j] << std::endl;
        j++;
    }

    gROOT->SetBatch(kTRUE);
    int n = e_field.size();

    double norm=y[0];
    for(int i=0;i<n;i++)
    {
        y[i]=y[i]/norm;
    }

    TCanvas* canvas = new TCanvas("canvas", "LY", 800, 600);
    TGraph* graph = new TGraph(n, e_field.data(), y.data());
    graph->SetTitle("Light Yield ; Efield [kv_cm] ; LY ");
   
    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end())); 
    graph->GetYaxis()->SetRangeUser(0.9*(*std::min_element(y.begin(), y.end())), 1.1*(*std::max_element(y.begin(), y.end())));
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(("./figures/" + std::to_string(ch) + "/ly_ch" + std::to_string(ch) + ".png").c_str());

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
}

void waveforms_Analyser::plot_ly_avg(std::vector<int> list_ch)
{
    std::vector<std::vector<double>> vector;
    for(const int& this_ch : list_ch)
    {
        vector.push_back(this->get_ly_vector(this_ch));
    }

    std::size_t n = vector.size(); //number of channels
    std::size_t dim = vector[0].size(); //nummber of voltage
    std::vector<double> mean(dim, 0.0);
    std::vector<double> stddev(dim, 0.0);

   
    for (const std::vector<double>& this_vector : vector) 
    {
        for (std::size_t i = 0; i < dim; ++i) 
        {
            mean[i] += this_vector[i];
        }
    }
    for (std::size_t i = 0; i < dim; ++i)
    {
        mean[i] /= n;
    }

    // Calcula o desvio padrão
    for (const std::vector<double>& this_vector : vector) 
    {
        for (std::size_t i = 0; i < dim; ++i)
        {
            double diff = this_vector[i] - mean[i];
            stddev[i] += diff * diff;
        }
    }
    for (std::size_t i = 0; i < dim; ++i)
    {
        stddev[i] = std::sqrt(stddev[i] / n); // para amostral: dividir por (n - 1)
    }
  
   
    auto aux_analyser = this->get_wf_Analyser_vector_by_ch(list_ch[0]);
    std::vector<double> e_field(dim);
    for(int i=0; i<dim ; i++)
    {
        e_field[i] = aux_analyser[i]->get_voltage();  
    }

    std::vector<double> ex(dim, 0.0);     

    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("canvas", "LY", 800, 600);
    TGraphErrors* graph = new TGraphErrors(dim, e_field.data(), mean.data(), ex.data(), stddev.data());
    graph->SetTitle("LY; Efield [kv_cm]; Light Yield");

    graph->GetXaxis()->SetLimits(0 -10 , 10 + *std::max_element(e_field.begin(), e_field.end()));  
    double y_min = 1e9, y_max = -1e9;
    for (int i = 0; i < n; ++i)
    {
        double y1 = mean[i] - stddev[i];
        double y2 = mean[i] + stddev[i];
        if (y1 < y_min) y_min = y1;
        if (y2 > y_max) y_max = y2;
    }
   //graph->GetYaxis()->SetRangeUser(y_min , y_max );


    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP"); 

    canvas->SetGrid();
    canvas->Update();
    canvas->SaveAs(std::string("./figures/avg/ly.png").c_str()); 

    delete graph;
    delete canvas;

    gROOT->SetBatch(kFALSE);
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
