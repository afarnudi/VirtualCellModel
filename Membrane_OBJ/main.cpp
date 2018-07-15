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
    Trajectory.open(membrane.output_file_neme.c_str());
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
    
    Membrane  membrane("membrane_2D_mesh");
    for (int i=0;i<200;i++)
    {
        trajectory(membrane, "membrane");
    }

    
    
    return 0;
}

