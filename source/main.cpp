/** @file doxygen_example.cpp
 @author Lastname:Firstname:A00123456:cscxxxxx
 @version Revision 1.1
 @brief Illustrates doxygen-style comments for documenting a C++
 program file and the functions in that file.
 @details If you want to add any further detailed description of
 what is in the file, then place it here (after the first statement)
 and it will appear in the detailed description section of the HTML
 output description for the file.
 @date Monday, September 19, 2011
 */

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
#include "Chromatin.h"
#include "Actin.h"

#include "General_functions.hpp"
#include "ECM.hpp"
#include "write_functions.hpp"
#include "interaction.hpp"
#include "maps.hpp"
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
    double Buffer_temperature;
    double Bussi_tau;
    double Actin_Membrane_Bond_Coefficient;
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
    vector<Actin> Actins;
    
    bool Include_Membrane  = false;
    bool Include_Chromatin = false;
    bool Include_ECM       = false;
    bool Include_Actin     = false;
    
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
            if (GenConst::Num_of_Membranes == GenConst::Num_of_Chromatins) {
                ///put a flag for chromatin inside membrane
                Chromatins[i].import_config(chromatin_config_list[i], Membranes[i].return_min_radius_after_relaxation());
            } else {
                Chromatins[i].import_config(chromatin_config_list[i]);
            }
            
            Chromatins[i].generate_report();
        }
    }
    
    if (GenConst::Num_of_Actins!=0) {
        Include_Actin = true;
        Actins.resize(GenConst::Num_of_Actins);
        for (int i=0; i<GenConst::Num_of_Actins; i++) {
            Actins[i].set_file_time(buffer);
            Actins[i].set_index(i);
            Actins[i].import_config(actin_config_list[i]);
//            Actins[i].generate_report();
        }
    }
    
    if (Include_Actin && Include_Membrane) {
        for (int i=0; i<GenConst::Num_of_Actins; i++) {
            Actin_Membrane_shared_Node_Identifier(Actins[i], Membranes[i]);
            
            if (Membranes[i].return_relax_with_actin_flag()) {
                Membranes[i].Relax();
            }
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
//        Chromatins[0].shift_velocity(0, 0.05, 0);
    }
    
    
    int progress=0;
    cout<<"\nBeginnig the MD\nProgress:\n";
    
    for(int MD_Step=0 ;MD_Step<=GenConst::MD_num_of_steps ; MD_Step++){
        
        
        //Velocity Verlet first stage
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
                Chromatins[i].Force_Calculator_2();
            }
        }
        if (Include_Actin)
        {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_beginning(GenConst::MD_Time_Step);
                Actins[i].Elastic_Force_Calculator();
            }
        }
        
        //Shared Forces
        
        if (Include_Chromatin && Include_Membrane) {
            if (MD_Step%2000==0) {
                for (int i=0; i<Chromatins.size(); i++) {
                    Chromatin_Membrane_neighbour_finder(Chromatins[i], Membranes[i]);
                    Chromatin_Membrane_hard_sphere(Chromatins[i], Membranes[i]);
                }
            }

            
        }
        
        if (Include_Membrane && Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                Actin_Membrane_shared_Node_Identifier(Actins[i], Membranes[i]);
            }
        }
        
        
        
        
        
        //Velocity Verlet second stage
        if (Include_Membrane) {
            for (int i=0; i<Membranes.size(); i++) {
                Membranes[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {//I think we can remove the first clause.
                    Membranes[i].Thermostat_Bussi(GenConst::Buffer_temperature);
                }
                
            }
        }
        if (Include_Chromatin) {
            for (int i=0; i<Chromatins.size(); i++) {
                Chromatins[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
                    Chromatins[i].Thermostat_Bussi(GenConst::MD_T*0.01);
                }
            }
        }
        if (Include_Actin) {
            for (int i=0; i<Actins.size(); i++) {
                Actins[i].MD_Evolution_end(GenConst::MD_Time_Step);
                if (GenConst::MD_thrmo_step!=0 && MD_Step%GenConst::MD_thrmo_step==0 && MD_Step>1000) {
                    Actins[i].Thermostat_Bussi(GenConst::MD_T);
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

