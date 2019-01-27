#include "Membrane.h"
#include "General_functions.hpp"
//#include "Bussi_Thermostat.hpp"

void Membrane::omega_calculator(){

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
        double len=vector_length_squared(a);
        COM_omega[0]+=temp_cross[0]/len;
        COM_omega[1]+=temp_cross[1]/len;
        COM_omega[2]+=temp_cross[2]/len;
    }
    Omega[0]=COM_omega[0];
    Omega[1]=COM_omega[1];
    Omega[2]=COM_omega[2];
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i][0]+=COM_velocity[0];//-temp_cross[0]/Num_of_Nodes;
        Node_Velocity[i][1]+=COM_velocity[1];//-temp_cross[1]/Num_of_Nodes;
        Node_Velocity[i][2]+=COM_velocity[2];//-temp_cross[2]/Num_of_Nodes;
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
    }
    
}

void Membrane::omega_calculator_2(){
    
    update_COM_velocity();
    update_COM_position();
    
    double COM_omega[3]={0,0,0};
    double COM_angular_momentum[3]={0,0,0};
    double temp_cross[3]={0};
    double Moment_of_inertia_COM[3][3]={0};
    double Moment_of_inertia_COM_inverse[3][3];
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
        
        Moment_of_inertia_COM[0][0]+=Node_Position[i][1]*Node_Position[i][1]+Node_Position[i][2]*Node_Position[i][2];
        Moment_of_inertia_COM[1][1]+=Node_Position[i][0]*Node_Position[i][0]+Node_Position[i][2]*Node_Position[i][2];
        Moment_of_inertia_COM[2][2]+=Node_Position[i][0]*Node_Position[i][0]+Node_Position[i][1]*Node_Position[i][1];
        Moment_of_inertia_COM[0][1]-=Node_Position[i][0]*Node_Position[i][1];
        Moment_of_inertia_COM[2][1]-=Node_Position[i][2]*Node_Position[i][1];
        Moment_of_inertia_COM[2][0]-=Node_Position[i][2]*Node_Position[i][0];
        Moment_of_inertia_COM[1][0]=Moment_of_inertia_COM[0][1];
        Moment_of_inertia_COM[1][2]=Moment_of_inertia_COM[2][1];
        Moment_of_inertia_COM[0][2]=Moment_of_inertia_COM[2][0];
        
        COM_angular_momentum[0]+=temp_cross[0];
        COM_angular_momentum[1]+=temp_cross[1];
        COM_angular_momentum[2]+=temp_cross[2];
    }
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            Moment_of_inertia_COM_inverse[i][j]=Moment_of_inertia_COM_inverse[i][j];
        }
    }
//    std::copy(&c[0][0], &Moment_of_inertia_COM[0][0]+9, &Moment_of_inertia_COM_inverse);
    matrix_inverse(Moment_of_inertia_COM_inverse);
    COM_omega[0]=Moment_of_inertia_COM_inverse[0][0]*COM_angular_momentum[0]+
                 Moment_of_inertia_COM_inverse[0][1]*COM_angular_momentum[1]+
                 Moment_of_inertia_COM_inverse[0][2]*COM_angular_momentum[2];
    COM_omega[1]=Moment_of_inertia_COM_inverse[1][0]*COM_angular_momentum[0]+
                 Moment_of_inertia_COM_inverse[1][1]*COM_angular_momentum[1]+
                 Moment_of_inertia_COM_inverse[1][2]*COM_angular_momentum[2];
    COM_omega[2]=Moment_of_inertia_COM_inverse[2][0]*COM_angular_momentum[0]+
                 Moment_of_inertia_COM_inverse[2][1]*COM_angular_momentum[1]+
                 Moment_of_inertia_COM_inverse[2][2]*COM_angular_momentum[2];
    
    double temp_k_omega=0.5*(COM_omega[0]*COM_angular_momentum[0] +
                             COM_omega[1]*COM_angular_momentum[1] +
                             COM_omega[2]*COM_angular_momentum[2]);
    
    delta_k_angular=temp_k_omega-k_angular;
    k_angular=temp_k_omega;
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i][0]+=COM_velocity[0];
        Node_Velocity[i][1]+=COM_velocity[1];
        Node_Velocity[i][2]+=COM_velocity[2];
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
    }
}
