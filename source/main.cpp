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
//    ECM surface("ECM_substrait_2", 0, -5, 0);
    ECM bacteria("ECM_Bacteria", 0, 4, 0);
    Membrane membrane("Membrane_Bacteria", 0, 5, 0);
    membrane.shift_position(0, 9.5, 0);
//    membrane.shift_velocity(0, -0.1, 0);
    //    Membrane  membrane("new_membrane");
    //    Membrane particle("newparticle");
    
    //begining of MD loop
    cout<<"Beginnig the MD\n";
    vector<int> bacteria_membrane_neighbour_list;
    bacteria_membrane_neighbour_list.resize(membrane.return_num_of_nodes(), -1);
//    vector<int> surface_membrane_neighbour_list;
//    surface_membrane_neighbour_list.resize(surface.return_num_of_nodes(), -1);
    bool costume_interaction_flag=false;
//    surface.set_interaction_range(2.0);
    bacteria.set_interaction_range(4.0);
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++){
        
        membrane.Membrane_MD_Evolution();
        membrane.Elastic_Force_Calculator();
//        interaction_1(MD_Step, membrane, surface, surface_membrane_neighbour_list);
        interaction_2(MD_Step, membrane, bacteria, bacteria_membrane_neighbour_list, costume_interaction_flag);
        interaction_3(membrane);
//        Node_ecm_Barrier(membrane, bacteria, bacteria_membrane_neighbour_list);
//        cout<<MD_Step<<endl;
//            membrane.ConstantSurfaceForceLocalTriangles();
//             if(MD_Step!=0 && MD_Step%100==0) // thermostate
//                {
//                 //cout<<"mem ";
//        //                Thermostat_2(membrane);
//                //cout<<"new ";
//        //        Thermostat_2(particle);
//                cout<<(MD_Step*100/MD_num_of_steps)<<"% done"<<endl;
//                }
        if (MD_Step%100==0) //saving Results for membrane
        {
//            membrane.calculate_average_force();
//            Results(bacteria, "bacteria", buffer);
//            Results(surface, "surface", buffer);
            Results(membrane, "membrane", buffer);
            double percent=100*MD_Step/MD_num_of_steps;
            cout<<percent<<endl;
            
        }
        //    particle.Membrane_MD_Evolution();
        //    particle.Elastic_Force_Calculator();
        //    particle.ConstantSurfaceForceLocalTriangles();
        
        //    if (MD_Step%200==0)//saving results for particle
        //        {
        //        Results(particle, "particle" ,buffer);
        //        } //End of if (MD_Step%200==0)
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
        V_com[0]+=membrane.Membrane_Node_Velocity[i][0];
        V_com[1]+=membrane.Membrane_Node_Velocity[i][1];
        V_com[2]+=membrane.Membrane_Node_Velocity[i][2];
    }
    
    
    V_com[0]/=membrane.return_num_of_nodes();
    V_com[2]/=membrane.return_num_of_nodes();
    V_com[1]/=membrane.return_num_of_nodes();
    
    //    cout<<"COM before= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\nV_X="<<V_com[0]<<"\tV_Y="<<V_com[1]<<"\tV_Z="<<V_com[2]<<endl;
    
    
    double alpha;
    //----------------------membrane---------------------
    
    //    alpha=sqrt(  (3*Membrane_num_of_Nodes*KT) / kineticenergymembrane( membrane.Membrane_Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( membrane.Membrane_Node_Velocity,vnuclei
    double Kinetic_energy=0;
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        membrane.Membrane_Node_Velocity[i][0]-=V_com[0];
        membrane.Membrane_Node_Velocity[i][1]-=V_com[1];
        membrane.Membrane_Node_Velocity[i][2]-=V_com[2];
    }
    
    for (int i=0; i<membrane.return_num_of_nodes(); i++) {
        Kinetic_energy+=membrane.Membrane_Node_Velocity[i][0]*membrane.Membrane_Node_Velocity[i][0]+membrane.Membrane_Node_Velocity[i][1]*membrane.Membrane_Node_Velocity[i][1]+membrane.Membrane_Node_Velocity[i][2]*membrane.Membrane_Node_Velocity[i][2];
    }
    Kinetic_energy*=membrane.Membrane_Node_Mass;
    
    
    
    alpha=sqrt(3*(membrane.return_num_of_nodes())*KT)/Kinetic_energy;
    //cout<<alpha<<endl;
    //cout<<V_com[0]<<"\t"<<V_com[1]<<"\t"<<V_com[2]<<"\n";
    
    for(int i=0;i<membrane.return_num_of_nodes();i++)
    {
        membrane.Membrane_Node_Velocity [i][0]*=alpha;
        membrane.Membrane_Node_Velocity [i][1]*=alpha;
        membrane.Membrane_Node_Velocity [i][2]*=alpha;
        membrane.Membrane_Node_Velocity[i][0]+=V_com[0];
        membrane.Membrane_Node_Velocity[i][1]+=V_com[1];
        membrane.Membrane_Node_Velocity[i][2]+=V_com[2];
    }
}
