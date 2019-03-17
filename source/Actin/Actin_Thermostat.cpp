#include "Actin.h"
#include "General_functions.hpp"
#include "Bussi_Thermostat.hpp"

void Actin::Thermostat_Bussi(double MD_T){
    update_COM_velocity();
    
    double alpha=0;
    Total_Kinetic_Energy=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Velocity[i][0]-=COM_velocity[0];
        Node_Velocity[i][1]-=COM_velocity[1];
        Node_Velocity[i][2]-=COM_velocity[2];
        Total_Kinetic_Energy+=Node_Velocity[i][0]*Node_Velocity[i][0]+Node_Velocity[i][1]*Node_Velocity[i][1]+Node_Velocity[i][2]*Node_Velocity[i][2];
    }
    double Num_degrees_of_freedom=3*Num_of_Nodes-3;
    Total_Kinetic_Energy*=0.5*Node_Mass;
    
    //    cout<<Total_Kinetic_Energy<<endl;
    //    T_Kinetic_Energy[MD_step%100]=Total_Kinetic_Energy;
    double sigma=0.5*Num_degrees_of_freedom*MD_T*GenConst::K;
    alpha=resamplekin(Total_Kinetic_Energy, sigma, Num_degrees_of_freedom, GenConst::Bussi_tau);
    alpha=sqrt(alpha/Total_Kinetic_Energy);
    //    cout<<"T = "<<2*Total_Kinetic_Energy/Num_degrees_of_freedom<<endl;
    //    cout<<"alpha_ Bussi = "<<alpha;
    //    cout<<endl;
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity [i][0]*=alpha;
        Node_Velocity [i][1]*=alpha;
        Node_Velocity [i][2]*=alpha;
        Node_Velocity [i][0]+=COM_velocity[0];
        Node_Velocity [i][1]+=COM_velocity[1];
        Node_Velocity [i][2]+=COM_velocity[2];
    }
}

