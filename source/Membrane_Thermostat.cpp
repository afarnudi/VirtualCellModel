#include "Membrane.h"
#include "General_functions.hpp"
#include "Bussi_Thermostat.hpp"



void Membrane::Thermostat_Bussi(double MD_T){
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
    alpha=resamplekin(Total_Kinetic_Energy, sigma, Num_degrees_of_freedom, GenConst::MD_thrmo_step);
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

void Membrane::Thermostat_2(double MD_T){
    update_COM_velocity();
    
    double alpha;
    Total_Kinetic_Energy=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Velocity[i][0]-=COM_velocity[0];
        Node_Velocity[i][1]-=COM_velocity[1];
        Node_Velocity[i][2]-=COM_velocity[2];
    
        Total_Kinetic_Energy+=Node_Velocity[i][0]*Node_Velocity[i][0]+Node_Velocity[i][1]*Node_Velocity[i][1]+Node_Velocity[i][2]*Node_Velocity[i][2];
    }
    Total_Kinetic_Energy*=Node_Mass/2;
    double Temperature=2*Total_Kinetic_Energy/(3*Num_of_Nodes-3);
    
    
    
    alpha=sqrt(GenConst::K*MD_T/Temperature);
    
//    cout<<"T = "<<Temperature<<endl;
//    cout<<"alpha_ thermo_2 = "<<alpha;
//    cout<<endl;
    
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

void Membrane::Thermostat_N6(double MD_T){
    
    update_COM_velocity();
    update_COM_position();
    
    double COM_omega[3]={0,0,0};
    double COM_angular_momentum[3]={0,0,0};
    double temp_cross[3]={0};
    double Moment_of_inertia_COM[3][3]={0};
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
        
        //
        Moment_of_inertia_COM[0][0]+=Node_Position[i][1]*Node_Position[i][1]+Node_Position[i][2]*Node_Position[i][2];
        Moment_of_inertia_COM[1][1]+=Node_Position[i][0]*Node_Position[i][0]+Node_Position[i][2]*Node_Position[i][2];
        Moment_of_inertia_COM[2][2]+=Node_Position[i][0]*Node_Position[i][0]+Node_Position[i][1]*Node_Position[i][1];
        Moment_of_inertia_COM[0][1]-=Node_Position[i][0]*Node_Position[i][1];
        Moment_of_inertia_COM[2][1]-=Node_Position[i][2]*Node_Position[i][1];
        Moment_of_inertia_COM[2][0]-=Node_Position[i][2]*Node_Position[i][0];
        Moment_of_inertia_COM[1][0]=Moment_of_inertia_COM[0][1];
        Moment_of_inertia_COM[1][2]=Moment_of_inertia_COM[2][1];
        Moment_of_inertia_COM[0][2]=Moment_of_inertia_COM[2][0];
        //
        
        COM_angular_momentum[0]+=temp_cross[0];
        COM_angular_momentum[1]+=temp_cross[1];
        COM_angular_momentum[2]+=temp_cross[2];
        
    }
    matrix_inverse(Moment_of_inertia_COM);
    COM_omega[0]=Moment_of_inertia_COM[0][0]*COM_angular_momentum[0]+Moment_of_inertia_COM[0][1]*COM_angular_momentum[1]+Moment_of_inertia_COM[0][2]*COM_angular_momentum[2];
    COM_omega[1]=Moment_of_inertia_COM[1][0]*COM_angular_momentum[0]+Moment_of_inertia_COM[1][1]*COM_angular_momentum[1]+Moment_of_inertia_COM[1][2]*COM_angular_momentum[2];
    COM_omega[2]=Moment_of_inertia_COM[2][0]*COM_angular_momentum[0]+Moment_of_inertia_COM[2][1]*COM_angular_momentum[1]+Moment_of_inertia_COM[2][2]*COM_angular_momentum[2];
    
    Total_Kinetic_Energy=0;
    
    
    for (int i=0; i<Num_of_Nodes; i++) {
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0] -= temp_cross[0];
        Node_Velocity[i][1] -= temp_cross[1];
        Node_Velocity[i][2] -= temp_cross[2];
        double b[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
        Total_Kinetic_Energy+=vector_length_squared(b);
    }
    double Num_degrees_of_freedom=3*Num_of_Nodes-6;
    double sigma=0.5*Num_degrees_of_freedom*MD_T*GenConst::K;
    double alpha=resamplekin(Total_Kinetic_Energy, sigma, Num_degrees_of_freedom, GenConst::MD_thrmo_step);
    alpha=sqrt(alpha/Total_Kinetic_Energy);
    
    Total_Kinetic_Energy*=0.5*Node_Mass;
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity [i][0]*=alpha;
        Node_Velocity [i][1]*=alpha;
        Node_Velocity [i][2]*=alpha;
        
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0]+=COM_velocity[0] + temp_cross[0];
        Node_Velocity[i][1]+=COM_velocity[1] + temp_cross[1];
        Node_Velocity[i][2]+=COM_velocity[2] + temp_cross[2];
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
        
        
    }
}
