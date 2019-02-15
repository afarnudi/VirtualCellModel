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
#include "Global_functions.hpp"

namespace GenConst {
    int MD_num_of_steps;
    int MD_traj_save_step;
    double MD_Time_Step;
    double MD_T;
    double K;
    int MD_thrmo_step;
    int MC_step;
    double Mem_fluidity;
    double Lbox;
    bool Periodic_condtion_status;
    int Num_of_Membranes;
    int Num_of_Chromatins;
    int Num_of_Actins;
    string trajectory_file_name;
    bool File_header;
    bool Relaxation;
    double Buffer_temperature;
    double Bussi_tau;
}


int main(int argc, char **argv)
{
    //time
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    string general_file_name="general-config.txt";
    vector<string> membrane_config_list;
    vector<string> chromatin_config_list;
    vector<string> actin_config_list;
    
    read_general_parameters(general_file_name, membrane_config_list, chromatin_config_list, actin_config_list);
    
    ofstream Trajectory;
    string traj_file_name="Results/"+GenConst::trajectory_file_name+buffer+".xyz";
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    
    vector<Membrane> Membranes;
    vector<Chromatin> Chromatins;
    
    bool Include_Membrane  = false;
    bool Include_Chromatin = false;
    bool Include_ECM       = false;
    
    if (GenConst::Num_of_Membranes!=0) {
        Include_Membrane = true;
        Membranes.resize(GenConst::Num_of_Membranes);
        for (int i=0; i<GenConst::Num_of_Membranes; i++) {
            Membranes[i].set_file_time(buffer);
            Membranes[i].set_index(i);
            Membranes[i].import_config(membrane_config_list[i]);
            Membranes[i].generate_report();
        }
    }
    
    if (GenConst::Num_of_Chromatins!=0) {
        Include_Chromatin = true;
        Chromatins.resize(GenConst::Num_of_Chromatins);
        for (int i=0; i<GenConst::Num_of_Chromatins; i++) {
            Chromatins[i].set_file_time(buffer);
            Chromatins[i].set_index(i);
            Chromatins[i].import_config(chromatin_config_list[i]);
            Chromatins[i].generate_report();
        }
    }
    
    
    
    if (Include_ECM)
    {
        ECM ecm("ECM",0,0,0);
        //        EcmGeneratingReport(buffer, ecm);
    }
    
    
    int num_of_elements=0;
    if (Include_Membrane) {
        for (int i=0; i<Membranes.size(); i++) {
            num_of_elements+=Membranes[i].return_num_of_nodes();
        }
    }
    if (Include_Chromatin) {
        for (int i=0; i<Chromatins.size(); i++) {
            num_of_elements+=Chromatins[i].return_num_of_nodes();
        }
        Chromatins[0].shift_velocity(0, 0.05, 0);
    }
    
    
    int progress=0;
    cout<<"\nBeginnig the MD\nProgress:\n";
    
    for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_steps ; MD_Step++){
        
        
        
        if (Include_Membrane)
        {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
                Membranes[i].Elastic_Force_Calculator(0);
            }
        }
        if (Include_Chromatin)
        {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
                Chromatins[i].Elastic_Force_Calculator();
            }
        }
        
        
        if (Include_Chromatin && Include_Membrane) {
            
            if (MD_Step%2000==0) {
                for (int i=0; i<Chromatins.size(); i++) {
                    for (int j=0; j<Membranes.size(); j++) {
                        Chromatin_Membrane_neighbour_finder(Chromatins[i], Membranes[j]);
                    }
                }
            }
            
            for (int i=0; i<Chromatins.size(); i++) {
                for (int j=0; j<Membranes.size(); j++) {
                    Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[j]);
                }
            }
            
        }
        
        if (Include_Membrane) {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {//I think we can remove the first clause.
//                    set_temperature(MD_Step, GenConst::MD_T, 1000000);
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
//                    Membranes[i].Thermostat_2(GenConst::MD_T);
//                    Membranes[i].Thermostat_N6(GenConst::MD_T);
                    Membranes[i].write_parameters(MD_Step);
                    
                }
                
            }
        }
        if (Include_Chromatin) {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.1);
                }
            }
        }
        
        
        
        if (MD_Step%GenConst::MD_traj_save_step==0) //saving Results
        {
            Trajectory << num_of_elements<<endl;
            Trajectory << " nodes  "<<endl;
            
            
            if (Include_Membrane) {
                for (int i=0; i<Membranes.size(); i++) {
                    string label="Membrane_"+to_string(i);
                    Membranes[i].write_traj(traj_file_name, label);
                    Membranes[i].export_for_resume(MD_Step);
                }
            }
            
            if (Include_Chromatin) {
                for (int i=0; i<Chromatins.size(); i++) {
                    string label="Chromatin_"+to_string(i);
                    Chromatins[i].write_traj(traj_file_name, label);
                    Chromatins[i].export_for_resume(MD_Step);
                }
            }
        }// End of if (MD_Step%100==0)
        
        
        if (int(100*MD_Step/GenConst::MD_num_of_steps)>progress){
            cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;
            progress+=5;
        }
        
        
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    cout<<"[ 100% ]\t step: "<<GenConst::MD_num_of_steps<<"\n";
    cout<<"\nDone!"<<endl;
    printf("Time taken: %.2f Minutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
    return 0;
}


