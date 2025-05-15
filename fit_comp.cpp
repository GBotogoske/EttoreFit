#include <iostream>
#include <string>
#include <vector>

#include "waveform_Analyser.hh"
#include "waveforms_Analyser.hh"

#include "utils.hh"

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

int main(int argc, char** argv)
{
    TApplication app("app", &argc, argv); // Necessário para exibir o canvas

    /* string file_path = "/home/gabriel/Documents/protodune/waveforms/wf_ch9_0kV.txt";
    string template_path = "template/ch_25_endpoint_112_avg.txt";

    waveform_Analyser* analyzer = new waveform_Analyser(file_path);
    analyzer->set_template(template_path);
    analyzer->set_ch_and_voltage_file();
    */
    //---------------------------------------------
    waveforms_Analyser* group_analyser = new waveforms_Analyser();
    group_analyser->set_folder("/home/gabriel/Documents/protodune/waveforms/");
    group_analyser->set_template_file("/home/gabriel/Documents/protodune/efield/template_list.txt");
    group_analyser->update();  

    //---------------------------------------------
 
    //vamos testar o FIT
    std::vector<double> par(6,0.0);
    par[0] = -1;
    par[1] = -0.1;
    par[2] = 100e-9;
    par[3] = 1500e-9;
    par[4] = 0;
    par[5] = -60;


    group_analyser->fit_all_channels(20,&par[0],6);


   std::vector<int> list_ch;
   list_ch.push_back(0);
   list_ch.push_back(1);
   list_ch.push_back(8);
   list_ch.push_back(9);
   list_ch.push_back(10);
   list_ch.push_back(11);
   list_ch.push_back(18);
   list_ch.push_back(19);
   

   group_analyser->plot_all_avg(list_ch);
    //group_analyser->plot_all_channels(20);

  /*   group_analyser->fit_channel(9,&par[0],6);
    group_analyser->plot_all(9);
 *//* 
    par[0] = -30;
    par[1] = -1;
    par[2] = 1e-9;
    par[3] = 1500e-9;
    par[4] = 0;
    par[5] = 55;
 */
    //group_analyser->fit_channel(9,&par[0],6);
    //group_analyser->plot_all(9);
    

    //analyzer->fit_0(&par[0],6);
    //std::cout << "number of times that calculated the convolution: " << analyzer->cont_conv << std::endl; 

    app.Run();

    delete group_analyser;

    return 0;
}