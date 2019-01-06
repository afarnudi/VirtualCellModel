#include "Chromatin.h"

void Chromatin::Thermostat_2(double MD_KT){
    double V_com[3];
    
    V_com[0]=0;
    V_com[1]=0;
    V_com[2]=0;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        V_com[0]+=Node_Velocity[i][0];
        V_com[1]+=Node_Velocity[i][1];
        V_com[2]+=Node_Velocity[i][2];
    }
    
    
    V_com[0]/=Num_of_Nodes;
    V_com[2]/=Num_of_Nodes;
    V_com[1]/=Num_of_Nodes;
    
    double alpha;
    
    //    alpha=sqrt(  (3*Membrane_num_of_Nodes*KT) / kineticenergymembrane( membrane.Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( membrane.Node_Velocity,vnuclei
    double Kinetic_energy=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Velocity[i][0]-=V_com[0];
        Node_Velocity[i][1]-=V_com[1];
        Node_Velocity[i][2]-=V_com[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        Kinetic_energy+=Node_Velocity[i][0]*Node_Velocity[i][0]+Node_Velocity[i][1]*Node_Velocity[i][1]+Node_Velocity[i][2]*Node_Velocity[i][2];
    }
    Kinetic_energy*=Node_Mass;
    
    
    
    alpha=sqrt(3*Num_of_Nodes*MD_KT/Kinetic_energy);
//    cout<<"3*Num_of_Nodes*MD_KT="<<3*Num_of_Nodes*MD_KT<<"\tKinetic_energy= "<<Kinetic_energy<<endl;
//    cout<<"V_com= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\talpha= "<<alpha<<endl;
//    cout<<alpha<<endl;
    //cout<<V_com[0]<<"\t"<<V_com[1]<<"\t"<<V_com[2]<<"\n";
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity [i][0]*=alpha;
        Node_Velocity [i][1]*=alpha;
        Node_Velocity [i][2]*=alpha;
        Node_Velocity[i][0]+=V_com[0];
        Node_Velocity[i][1]+=V_com[1];
        Node_Velocity[i][2]+=V_com[2];
    }
}
