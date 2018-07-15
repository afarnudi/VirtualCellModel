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



void Results (Membrane membrane, string label, char* buffer)
{
	//output file names:
	
	string energy_file_name;
	string traj_file_name;
	membrane.output_file_neme+=buffer;
	traj_file_name="Results_";
	traj_file_name+=membrane.output_file_neme;
	traj_file_name+=".xyz";
	energy_file_name="Results_Potential_Energy";
	energy_file_name+=membrane.output_file_neme;
	
	//trajectory:
	
    ofstream Trajectory;

    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory << membrane.Membrane_num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< membrane.Membrane_num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<membrane.Membrane_Node_Position[j][0]<< setw(20)<<membrane.Membrane_Node_Position[j][1]<< setw(20)<<membrane.Membrane_Node_Position[j][2]<<endl;
    }
	
	//Energy:
	ofstream Potential_Energy;
	Potential_Energy.open(energy_file_name.c_str(), ios::app);
	Potential_Energy<<membrane.Total_Potential_Energy<<endl;
}

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
	Membrane  membrane("new_membrane");
	Membrane particle("newparticle");  

	//begining of MD loop
	for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    {
    membrane.Membrane_MD_Evolution();
	membrane.Elastic_Force_Calculator();
	membrane.ConstantSurfaceForceLocalTriangles();
	 if(MD_Step!=0 && MD_Step%100==0) // thermostate
        {
			//cout<<"mem ";
            Thermostat_2(membrane);
			//cout<<"new ";
			Thermostat_2(particle);
        }
	if (MD_Step%200==0) //saving Results for membrane
		{
		Results(membrane, "membrane", buffer);
		}
    particle.Membrane_MD_Evolution();
	particle.Elastic_Force_Calculator();
	particle.ConstantSurfaceForceLocalTriangles();
	
	if (MD_Step%200==0)//saving results for particle 
		{
		Results(particle, "particle" ,buffer);
		} //End of if (MD_Step%200==0)
	} //End of for (int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
	/*
	
	
    //cheking the Membrane_triangle_pair_and_Edges_identifier
	cout<<membrane.Membrane_num_of_Triangle_Pairs<<endl;
	ofstream Edges;
    Edges.open("Edges");
	for (int i=0; i<membrane.Membrane_num_of_Triangle_Pairs ; i++)
	{
		 Edges << "TNP" <<setprecision(5)<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][0]<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][1]<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][2]<<setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][3]<<setw(10)<<"ME"<<setw(10)<<membrane.Membrane_Edges[i][0]<<setw(10)<<membrane.Membrane_Edges[i][1]<<endl;
	}
    */

    return 0;
}

void Thermostat_2(Membrane membrane)
{
    double V_com[3];
    
    V_com[0]=0;
    V_com[1]=0;
    V_com[2]=0;
    for (int i=0; i<membrane.Membrane_num_of_Nodes; i++) {
        V_com[0]+=membrane.Membrane_Node_Velocity[i][0];
        V_com[1]+=membrane.Membrane_Node_Velocity[i][1];
        V_com[2]+=membrane.Membrane_Node_Velocity[i][2];
    }
    
    
    V_com[0]/=membrane.Membrane_num_of_Nodes;
    V_com[2]/=membrane.Membrane_num_of_Nodes;
    V_com[1]/=membrane.Membrane_num_of_Nodes;
    
//    cout<<"COM before= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\nV_X="<<V_com[0]<<"\tV_Y="<<V_com[1]<<"\tV_Z="<<V_com[2]<<endl;
    
    
    double alpha;
    //----------------------membrane---------------------
    
    //    alpha=sqrt(  (3*Membrane_num_of_Nodes*KT) / kineticenergymembrane( membrane.Membrane_Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( membrane.Membrane_Node_Velocity,vnuclei
    double Kinetic_energy=0;
    for (int i=0; i<membrane.Membrane_num_of_Nodes; i++) {
        membrane.Membrane_Node_Velocity[i][0]-=V_com[0];
        membrane.Membrane_Node_Velocity[i][1]-=V_com[1];
        membrane.Membrane_Node_Velocity[i][2]-=V_com[2];
    }
    
    for (int i=0; i<membrane.Membrane_num_of_Nodes; i++) {
        Kinetic_energy+=membrane.Membrane_Node_Velocity[i][0]*membrane.Membrane_Node_Velocity[i][0]+membrane.Membrane_Node_Velocity[i][1]*membrane.Membrane_Node_Velocity[i][1]+membrane.Membrane_Node_Velocity[i][2]*membrane.Membrane_Node_Velocity[i][2];
    }
    Kinetic_energy*=membrane.Membrane_Node_Mass;
    


    alpha=sqrt(3*(membrane.Membrane_num_of_Nodes)*KT)/Kinetic_energy;
	//cout<<alpha<<endl;
	//cout<<V_com[0]<<"\t"<<V_com[1]<<"\t"<<V_com[2]<<"\n";
    
    for(int i=0;i<membrane.Membrane_num_of_Nodes;i++)
    {
        membrane.Membrane_Node_Velocity [i][0]*=alpha;
        membrane.Membrane_Node_Velocity [i][1]*=alpha;
        membrane.Membrane_Node_Velocity [i][2]*=alpha;
        membrane.Membrane_Node_Velocity[i][0]+=V_com[0];
        membrane.Membrane_Node_Velocity[i][1]+=V_com[1];
        membrane.Membrane_Node_Velocity[i][2]+=V_com[2];
    }
}
