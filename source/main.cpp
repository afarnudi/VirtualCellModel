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

int MD_num_of_steps=300000;//35000// number of MD stps
int MD_traj_save_step=2000;//The step on which the trajector of the membrane is saved.
double MD_Time_Step=0.001; // time length of steps in MD
double MD_KT=1.0;  // KT the quanta of energy
int MD_thrmo_step=100; //
int MC_step=1;
double Mem_fluidity=0.002; //Used in the MC step
double Lbox=1000.0;///    (size of square periodic box-1)
bool Periodic_condtion_status=false; //status 0.0 = false (The Periodic update will not be executed in the 'Main MD' loop). status = 1.0 = true

void Thermostat_2(Membrane membrane, double MD_KT);
void read_general_parameters(map<string, double> general_param_map, string input_file_name);

int main(int argc, char **argv)
{
    //time
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_time_%H_%M",now);
    
    
    //indicating what kind of objects participate in simulation (we can set the value this flags via an input file)
    map<string, double> general_param_map;
    string general_file_name="test_map.txt";
    read_general_parameters(general_param_map, general_file_name);
    
    bool Include_Membrane =true;
    bool Include_ECM= false;
    bool Include_Particle=false;
    bool resume=false;
    //initialling Membrane classes
    //if (Include_Membrane)
    //this does not work because the varibles of type membrane and ecm are defined just in if statement scope. we have to find a solution for this if we want to use input file.
    Membrane membrane;
    if (resume) {
        membrane.import("Resume_2018_10_14_21_31.txt");
    } else{
        membrane.initialise("Membrane", 0, 0, 0);
    }
    
    cout<<membrane.return_num_of_nodes()<<endl;
    
    if (Include_Membrane)
    {
        MembraneGeneratingReport(buffer,membrane);
    }
    if (Include_ECM)
    {
        ECM ecm("ECM",0,0,0);
        EcmGeneratingReport(buffer, ecm);
    }
    if (Include_Particle)
    {
        Membrane Particle;
        Particle.initialise("particle",0,0,0);
        ParticleGeneratingReport(buffer, Particle);
    }
    
    //    }
    //begining of MD loop
    cout<<"Beginnig the MD\nProgress:\n";
    int progress=0;
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++){
        
        if (Include_Membrane)
        {
            membrane.MD_Evolution_beginning(MD_Time_Step);
            membrane.Elastic_Force_Calculator(0);
            membrane.MD_Evolution_end(MD_Time_Step);
            
            
            if (MD_Step%MD_traj_save_step==0) //saving Results for membrane
            {
                membrane.export_for_resume(buffer, MD_Step);
            }// End of if (MD_Step%100==0)
            if (MD_Step%100==0) //saving Results for membrane
            {
                Results(membrane, "membrane", buffer);
                //double percent=100*MD_Step/MD_num_of_steps;
                //cout<<percent<<endl;
            }// End of if (MD_Step%100==0)
            if (int(100*MD_Step/MD_num_of_steps)>progress){
                cout<<"[ "<<progress<<"% ]\t step: "<<MD_Step<<"\r" << std::flush;;
                progress+=5;
            }
            
        }//End of if (Include_Membrane==true)
        
    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    
    cout<<"done"<<endl;
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

void read_general_parameters(map<string, double> general_param_map, string input_file_name){
    ifstream read_map(input_file_name.c_str());
    string param_name;
    double param_value;
    //    map<string, double> general_param_map;
    map<string, double>::iterator it;
    
    while (read_map>>param_name>>param_value) {
        general_param_map[param_name]=param_value;
    }
    
    param_name="MD_num_of_steps";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MD_num_of_steps=it->second;
    }
    param_name="MD_traj_save_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MD_traj_save_step=it->second;
    }
    param_name="MD_Time_Step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MD_Time_Step=it->second;
    }
    param_name="MD_KT";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MD_KT=it->second;
    }
    param_name="MD_thrmo_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MD_thrmo_step=it->second;
    }
    param_name="MC_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        MC_step=it->second;
    }
    param_name="Mem_fluidity";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        Mem_fluidity=it->second;
    }
    param_name="Lbox";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        Lbox=it->second;
    }
    param_name="Periodic_condtion_status";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        if (it->second==0.0) {
            Periodic_condtion_status=false;
        } else {
            Periodic_condtion_status=true;
        }
    }
}
