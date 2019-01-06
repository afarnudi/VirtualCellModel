#include "Membrane.h"
#include "General_functions.hpp"

void Membrane::Thermostat_2(double MD_KT){
    update_COM_velocity();
    
    double alpha;
    double Kinetic_energy=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Velocity[i][0]-=COM_velocity[0];
        Node_Velocity[i][1]-=COM_velocity[1];
        Node_Velocity[i][2]-=COM_velocity[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        Kinetic_energy+=Node_Velocity[i][0]*Node_Velocity[i][0]+Node_Velocity[i][1]*Node_Velocity[i][1]+Node_Velocity[i][2]*Node_Velocity[i][2];
    }
    Kinetic_energy*=Node_Mass/(3*Num_of_Nodes);
    
    
    
    alpha=sqrt(MD_KT/Kinetic_energy);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity [i][0]*=alpha;
        Node_Velocity [i][1]*=alpha;
        Node_Velocity [i][2]*=alpha;
        Node_Velocity[i][0]+=COM_velocity[0];
        Node_Velocity[i][1]+=COM_velocity[1];
        Node_Velocity[i][2]+=COM_velocity[2];
    }
}

void Membrane::Thermostat_N6(double MD_KT){
    
    update_COM_velocity();
    update_COM_position();
    double COM_omega[3]={0};
    double temp_cross[3]={0};
    
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Velocity[i][0] -= COM_velocity[0];
        Node_Velocity[i][1] -= COM_velocity[1];
        Node_Velocity[i][2] -= COM_velocity[2];
        
        Node_Position[i][0] -= COM_position[0];
        Node_Position[i][1] -= COM_position[1];
        Node_Position[i][2] -= COM_position[2];
        
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        double b[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
        crossvector(temp_cross, a, b);
        double len=vector_length_squared(temp_cross);
        COM_omega[0]+=temp_cross[0]/len;
        COM_omega[1]+=temp_cross[1]/len;
        COM_omega[2]+=temp_cross[2]/len;
    }
    
    double Kinetic_energy=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        crossvector(temp_cross, COM_omega, a);
        Node_Velocity[i][0]-=temp_cross[0];
        Node_Velocity[i][1]-=temp_cross[1];
        Node_Velocity[i][2]-=temp_cross[2];
        double b[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
        Kinetic_energy+=vector_length_squared(b);
    }
    double initial_temperature=Kinetic_energy*Node_Mass/(3*Num_of_Nodes-6);
    
    double alpha=sqrt(MD_KT/initial_temperature);
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity [i][0]*=alpha;
        Node_Velocity [i][1]*=alpha;
        Node_Velocity [i][2]*=alpha;
        
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0]+=COM_velocity[0]+temp_cross[0];
        Node_Velocity[i][1]+=COM_velocity[1]+temp_cross[1];
        Node_Velocity[i][2]+=COM_velocity[2]+temp_cross[2];
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
        
        
    }
}
