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
}


void Thermostat_2(Membrane membrane, double MD_KT);


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
    read_general_parameters(general_file_name, membrane_config_list);
    vector<Membrane> Membranes;
    bool Include_Membrane = false;
    if (GenConst::Num_of_Membranes!=0) {
        Include_Membrane = true;
        Membranes.resize(GenConst::Num_of_Membranes);
        for (int i=0; i<GenConst::Num_of_Membranes; i++) {
            Membranes[i].import_config(membrane_config_list[i]);
        }
    }
    
    bool Include_ECM= false;
    //    bool Include_Particle=false;
    bool resume=false;
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

void Thermostat_2(Membrane membrane, double MD_KT)
{
    double V_com[3];
    
    V_com[0]=0;
    V_com[1]=0;
    V_com[2]=0;
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        V_com[0]+=membrane.Node_Velocity[i][0];
        V_com[1]+=membrane.Node_Velocity[i][1];
        V_com[2]+=membrane.Node_Velocity[i][2];
    }
    
    
    V_com[0]/=membrane.return_num_of_nodes();
    V_com[2]/=membrane.return_num_of_nodes();
    V_com[1]/=membrane.return_num_of_nodes();
    
    //    cout<<"COM before= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\nV_X="<<V_com[0]<<"\tV_Y="<<V_com[1]<<"\tV_Z="<<V_com[2]<<endl;
    
    
    double alpha;
    //----------------------membrane---------------------
    
    //    alpha=sqrt(  (3*Membrane_num_of_Nodes*KT) / kineticenergymembrane( membrane.Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( membrane.Node_Velocity,vnuclei
    double Kinetic_energy=0;
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        membrane.Node_Velocity[i][0]-=V_com[0];
        membrane.Node_Velocity[i][1]-=V_com[1];
        membrane.Node_Velocity[i][2]-=V_com[2];
    }
    
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        Kinetic_energy+=membrane.Node_Velocity[i][0]*membrane.Node_Velocity[i][0]+membrane.Node_Velocity[i][1]*membrane.Node_Velocity[i][1]+membrane.Node_Velocity[i][2]*membrane.Node_Velocity[i][2];
    }
    Kinetic_energy*=membrane.Node_Mass;
    
    
    
    alpha=sqrt(3*(membrane.return_num_of_nodes())*MD_KT)/Kinetic_energy;
    //cout<<alpha<<endl;
    //cout<<V_com[0]<<"\t"<<V_com[1]<<"\t"<<V_com[2]<<"\n";
    
    for(int i=0;i<membrane.return_num_of_nodes();i++)
    {
        membrane.Node_Velocity [i][0]*=alpha;
        membrane.Node_Velocity [i][1]*=alpha;
        membrane.Node_Velocity [i][2]*=alpha;
        membrane.Node_Velocity[i][0]+=V_com[0];
        membrane.Node_Velocity[i][1]+=V_com[1];
        membrane.Node_Velocity[i][2]+=V_com[2];
    }
}


