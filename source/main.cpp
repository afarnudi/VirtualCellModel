#include <stdio.h>
#include <ctime>
#include <sstream>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>

#include "Membrane.h"
#include "General_functions.hpp"
#include "ECM.hpp"
#include "write_functions.hpp"
#include "interaction.hpp"
#include "maps.hpp"
#include "Chromatin.h"

namespace GenConst {
    int MD_num_of_steps;
    int MD_traj_save_step;
    double MD_Time_Step;
    double MD_KT;
    int MD_thrmo_step;
    int MC_step;
    double Mem_fluidity;
    double Lbox;
    bool Periodic_condtion_status;
    int Num_of_Membranes;
    int Num_of_Chromatins;
}


int main(int argc, char **argv)
{
    //time
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    
    //indicating what kind of objects participate in simulation (we can set the value this flags via an input file)
    //    map<string, double> general_param_map;
    string general_file_name="general-config.txt";
    vector<string> membrane_config_list;
    vector<string> chromatin_config_list;
    read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list);
    
    vector<Membrane> Membranes;
    vector<Chromatin> Chromatins;
    bool Include_Membrane = false;
    bool Include_Chromatin = false;
    
    if (GenConst::Num_of_Membranes!=0) {
        Include_Membrane = true;
        Membranes.resize(GenConst::Num_of_Membranes);
        for (int i=0; i<GenConst::Num_of_Membranes; i++) {
            Membranes[i].import_config(membrane_config_list[i]);
        }
    }
    
    if (GenConst::Num_of_Chromatins!=0) {
        Include_Chromatin = true;
        Chromatins.resize(GenConst::Num_of_Chromatins);
        for (int i=0; i<GenConst::Num_of_Chromatins; i++) {
            Chromatins[i].import_config(chromatin_config_list[i]);
        }
    }
    
    bool Include_ECM= false;
    //    bool Include_Particle=false;
    //    bool resume=false;
    //initialling Membrane classes
    //if (Include_Membrane)
    //this does not work because the varibles of type membrane and ecm are defined just in if statement scope. we have to find a solution for this if we want to use input file.
    //    Membrane membrane;
    //    if (resume) {
    //        membrane.import("Resume_2018_10_14_21_31.txt");
    //    } else{
    //        membrane.initialise("Membrane", 0, 0, 0);
    //    }
    for (int i=0; i<Membranes.size(); i++) {
        cout<<"Membrane_"<<i<<" # of nodes = "<<Membranes[i].return_num_of_nodes()<<endl;
    }
    
    
    if (Include_Membrane)
    {
        for (int i=0; i<Membranes.size(); i++) {
            Membranes[i].generate_report(buffer);
        }
        
    }
    if (Include_ECM)
    {
        ECM ecm("ECM",0,0,0);
        EcmGeneratingReport(buffer, ecm);
    }
    //    if (Include_Particle)
    //    {
    //        Membrane Particle;
    //        Particle.initialise("particle",0,0,0);
    //        ParticleGeneratingReport(buffer, Particle);
    //    }
    
    //    }
    //begining of MD loop
    cout<<"Beginnig the MD\nProgress:\n";
    int progress=0;
    for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_steps ; MD_Step++){
        
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
                Membranes[i].Elastic_Force_Calculator(0);
            }
            
            
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
                    Membranes[i].Thermostat_2(GenConst::MD_KT);
                }
            }
            
            
            
            if (MD_Step%GenConst::MD_traj_save_step==0) //saving Results for membrane
            {
                for (int i=0; i<Membranes.size(); i++) {
                    Membranes[i].export_for_resume(buffer, MD_Step);
                }
            }// End of if (MD_Step%100==0)
            if (MD_Step%100==0) //saving Results for membrane
            {
                for (int i=0; i<Membranes.size(); i++) {
                    string file_prefix="membrane_"+to_string(i);
                    Results(Membranes[i], file_prefix, buffer);
                }
                
                //double percent=100*MD_Step/MD_num_of_steps;
                //cout<<percent<<endl;
            }// End of if (MD_Step%100==0)
            if (int(100*MD_Step/GenConst::MD_num_of_steps)>progress){
                cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
                progress+=5;
            }
            
        }//End of if (Include_Membrane==true)
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    cout<<"[ 100% ]\t step: "<<GenConst::MD_num_of_steps<<"\n";
    cout<<"\nDone!"<<endl;
    printf("Time taken: %.2fminutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
    return 0;
}


