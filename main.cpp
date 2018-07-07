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



void trajectory (Membrane membrane, string label)
{
    ofstream Trajectory;
    Trajectory.open(membrane.output_file_neme.c_str(), ios::app);
    Trajectory << std:: fixed;
    Trajectory << membrane.Membrane_num_of_Nodes<<endl;
    Trajectory << " nodes  "<<endl;
    for(int j=0; j< membrane.Membrane_num_of_Nodes;j++) // saving trajectory
    {
        Trajectory << label <<setprecision(5)<< setw(20)<<membrane.Membrane_Node_Position[j][0]<< setw(20)<<membrane.Membrane_Node_Position[j][1]<< setw(20)<<membrane.Membrane_Node_Position[j][2]<<endl;
    }
}


int main(int argc, char **argv)
{
    /*clock_t tStart = clock();//Time the programme
     time_t t = time(0);   // get time now
     struct tm * now = localtime( & t );
     char buffer [80];
     strftime (buffer,80,"%Y_%m_%d_%H:%M",now);
     //outputfiles:
     /*string traj_file_name;
     
     traj_file_name="results/membrane_";
     traj_file_name +=buffer;
     traj_file_name +=".xyz";*/
    

    Membrane  membrane("new_membrane");
	Membrane particle("newparticle");
	for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    {
    membrane.Membrane_MD_Evolution();
	membrane.Elastic_Force_Calculator();
	trajectory(membrane, "membrane");
    particle.Membrane_MD_Evolution();
	particle.Elastic_Force_Calculator();
	trajectory(particle, "particle");	
	}
	/*
	for (int i=0;i<1; i++)
	{
	trajectory(membrane, "membrane");
	trajectory(particle, "particle");
    }
	*/
	
    //cheking the Membrane_triangle_pair_and_Edges_identifier
	cout<<membrane.Membrane_num_of_Triangle_Pairs<<endl;
	ofstream Edges;
    Edges.open("Edges");
	for (int i=0; i<membrane.Membrane_num_of_Triangle_Pairs ; i++)
	{
		 Edges << "TNP" <<setprecision(5)<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][0]<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][1]<< setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][2]<<setw(10)<<membrane.Membrane_Triangle_Pair_Nodes[i][3]<<setw(10)<<"ME"<<setw(10)<<membrane.Membrane_Edges[i][0]<<setw(10)<<membrane.Membrane_Edges[i][1]<<endl;
	}
    
    
    return 0;
}

