#include "Membrane.h"
#include "General_functions.hpp"
#include "Bussi_Thermostat.hpp"

void matrix_inverse (double mat[3][3]);

void Membrane::Thermostat_Bussi(double MD_KT){
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
    double sigma=0.5*Num_degrees_of_freedom*GenConst::MD_T;
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

void Membrane::Thermostat_2(double MD_KT){
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
    
    
    
    alpha=sqrt(MD_KT/Temperature);
    
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

void Membrane::Thermostat_N6(double MD_KT){
    
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
    
//    double initial_Omega[3]={COM_omega[0], COM_omega[1], COM_omega[2]};
//    double final_Omega[3]={0,0,0};
//    double final_ang_mom[3]={0,0,0};
//    double temp_I=0;
    Total_Kinetic_Energy=0;
    
    
    for (int i=0; i<Num_of_Nodes; i++) {
        double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0] -= temp_cross[0];
        Node_Velocity[i][1] -= temp_cross[1];
        Node_Velocity[i][2] -= temp_cross[2];
        double b[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
        Total_Kinetic_Energy+=vector_length_squared(b);
        
        
//        double cross_delete[3];
//        double c[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
//        crossvector(cross_delete, a, b);
//        temp_I+=vector_length_squared(a);
//        final_ang_mom[0]+=cross_delete[0];
//        final_ang_mom[1]+=cross_delete[1];
//        final_ang_mom[2]+=cross_delete[2];
    }
    
//    double alpha=sqrt(MD_KT/initial_temperature);
//    final_Omega[0]=Moment_of_inertia_COM[0][0]*final_ang_mom[0]+Moment_of_inertia_COM[0][1]*final_ang_mom[1]+Moment_of_inertia_COM[0][2]*final_ang_mom[2];
//    final_Omega[1]=Moment_of_inertia_COM[1][0]*final_ang_mom[0]+Moment_of_inertia_COM[1][1]*final_ang_mom[1]+Moment_of_inertia_COM[1][2]*final_ang_mom[2];
//    final_Omega[2]=Moment_of_inertia_COM[2][0]*final_ang_mom[0]+Moment_of_inertia_COM[2][1]*final_ang_mom[1]+Moment_of_inertia_COM[2][2]*final_ang_mom[2];
    
//    cout<<"initial Omega = "<<initial_Omega[0]<<"\t"<<initial_Omega[1]<<"\t"<<initial_Omega[2]<<"\n";
//    cout<<"  final Omega = "<<final_Omega[0]<<"\t"<<final_Omega[1]<<"\t"<<final_Omega[2]<<"\n";
//    cout<<" COM Velocity = "<<COM_velocity[0]<<"\t"<<COM_velocity[1]<<"\t"<<COM_velocity[2];
//    cout<<" COM Position = "<<COM_position[0]<<"\t"<<COM_position[1]<<"\t"<<COM_position[2];
//    cout<<endl;
    
    
    double Num_degrees_of_freedom=3*Num_of_Nodes-6;
    double sigma=0.5*Num_degrees_of_freedom*GenConst::MD_T;
    double alpha=resamplekin(Total_Kinetic_Energy, sigma, Num_degrees_of_freedom, GenConst::MD_thrmo_step);
    alpha=sqrt(alpha/Total_Kinetic_Energy);
//    double initial_temperature=Total_Kinetic_Energy*Node_Mass/(3*Num_of_Nodes-6);
    Total_Kinetic_Energy*=0.5*Node_Mass;
    
//    double alpha=sqrt(MD_KT/initial_temperature);
//    cout<<"alpha = "<<alpha<<endl;
    
    
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

void Membrane::omega(int MD_Step, double step){
//    string energy_file_name;
    string traj_file_name;
    
    traj_file_name="Results/Omega_"+GenConst::trajectory_file_name+"Membrane_"+to_string(mem_index)+"_"+file_time+".txt";
    //trajectory:
    
    ofstream Trajectory;
    
    Trajectory.open(traj_file_name.c_str(), ios::app);
    Trajectory << std:: fixed;
//    Trajectory <<Num_of_Nodes<<endl;
//    Trajectory << " nodes  "<<endl;
    
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
    Omega[0]+=COM_omega[0]/Num_of_Nodes;
    Omega[1]+=COM_omega[1]/Num_of_Nodes;
    Omega[2]+=COM_omega[2]/Num_of_Nodes;
    
    Trajectory<<MD_Step<<setw(20)<<COM_omega[0]<<setw(20)<<COM_omega[1]<<setw(20)<<COM_omega[2]<< setw(20)<<sqrt(COM_omega[0]*COM_omega[0]+COM_omega[2]*COM_omega[2]+COM_omega[1]*COM_omega[1])
    
    <<setw(20)<<step*Omega[0]/MD_Step<<setw(20)<<step*Omega[1]/MD_Step<<setw(20)<<step*Omega[2]/MD_Step<< setw(20)<<sqrt(Omega[0]*Omega[0]+Omega[2]*Omega[2]+Omega[1]*Omega[1])*step/MD_Step<<endl;

    for(int i=0;i<Num_of_Nodes;i++)
    {
        //            double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        //            double temp_cross[3]={0};
        //            crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0]+=COM_velocity[0];//-temp_cross[0]/Num_of_Nodes;
        Node_Velocity[i][1]+=COM_velocity[1];//-temp_cross[1]/Num_of_Nodes;
        Node_Velocity[i][2]+=COM_velocity[2];//-temp_cross[2]/Num_of_Nodes;
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
    }
    
}

void Membrane::equilibrate (void){
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
    Omega[0]+=COM_omega[0];
    Omega[1]+=COM_omega[1];
    Omega[2]+=COM_omega[2];
    
//    Trajectory<<MD_Step<< setw(20)<<step*Omega[0]/MD_Step<< setw(20)<<step*Omega[1]/MD_Step<< setw(20)<<step*Omega[2]/MD_Step<< setw(20)<<sqrt(Omega[0]*Omega[0]+Omega[2]*Omega[2]+Omega[1]*Omega[1])*step/MD_Step<<endl;
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        //            double a[3]={Node_Position[i][0], Node_Position[i][1], Node_Position[i][2]};
        //            double temp_cross[3]={0};
        //            crossvector(temp_cross, COM_omega, a);
        
        Node_Velocity[i][0]+=COM_velocity[0]-temp_cross[0]/Num_of_Nodes;
        Node_Velocity[i][1]+=COM_velocity[1]-temp_cross[1]/Num_of_Nodes;
        Node_Velocity[i][2]+=COM_velocity[2]-temp_cross[2]/Num_of_Nodes;
        
        Node_Position[i][0]+=COM_position[0];
        Node_Position[i][1]+=COM_position[1];
        Node_Position[i][2]+=COM_position[2];
    }
}

void matrix_inverse (double mat[3][3]){
    int i, j;
    float determinant = 0;
    
    //finding determinant
    for(i = 0; i < 3; i++)
        determinant += (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
    
    double mat_temp[3][3]={0};
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            mat_temp[i][j]=((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant;
        
    }
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            mat[i][j]=mat_temp[i][j];
        
    }
}
