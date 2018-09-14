#include <stdio.h>
#include "Membrane.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>
#include "General_functions.hpp"
#include "ECM.hpp"
#include "write_functions.hpp"
#include "interaction.hpp"



void Thermostat_2(Membrane membrane);

int main(int argc, char **argv)
{
    //time
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_%H:%M",now);
	
    //initialling Membrane classes
    Membrane membrane("new_membrane", 0, 0, 0);

    //begining of MD loop
    cout<<"Beginnig the MD\n";
	
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++){
        
        membrane.MD_Evolution();
        membrane.Elastic_Force_Calculator();
		
        if (MD_Step%100==0) //saving Results for membrane
        {
            Results(membrane, "membrane", buffer);
            double percent=100*MD_Step/MD_num_of_steps;
            cout<<percent<<endl;
            
        }

    } //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)

    cout<<"done"<<endl;
    printf("Time taken: %.2fminutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
    return 0;
}

void Thermostat_2(Membrane membrane)
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
    
    
    
    alpha=sqrt(3*(membrane.return_num_of_nodes())*KT)/Kinetic_energy;
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