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
    group_analyser->set_folder("./waveforms/");
    group_analyser->set_template_file("./template_list.txt");
    group_analyser->update();  

    //---------------------------------------------
 
    //vamos testar o FIT
    std::vector<double> par(6,0.0);
    par[0] = -3;
    par[1] = -0.1;
    par[2] = 50e-9;
    par[3] = 1600e-9;
    par[4] = 0;
    par[5] = -57;

    group_analyser->fit_all_channels(20,&par[0],6);

   std::vector<int> list_ch;
   for(int i=0;i<20;i++)
   {
      list_ch.push_back(i);   
   }
   
   group_analyser->plot_all_avg(list_ch);

    app.Run();

    delete group_analyser;

    return 0;
}