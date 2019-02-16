//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"
//void build_random_chain(void);

void Chromatin::initialise(void){
    cout<<"\nInitialising the Chromatin Class..."<<endl;
    build_random_chain();
    cout<<"Chromatin class initiated.\n";

}

void Chromatin::build_random_chain(void){
    
    double position_COM[3]={0};
    double velocity_COM[3]={0};
    
    Node_Position.resize(Num_of_Nodes);
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    
    
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
        Node_Position[i].resize(3,0);
    }
    const double PI  =3.141592653589793238463;
    
    for (int i=1; i<Num_of_Nodes; i++) {
        double theta=((double)rand()/(double)RAND_MAX)*PI;
        double phi=((double)rand()/(double)RAND_MAX)*2*PI;
//        double temp_x=Node_Position[i-1][0]+Node_radius*(2.1+((double)rand()/(double)RAND_MAX)*0.2)*sin(theta)*cos(phi);
//        double temp_y=Node_Position[i-1][1]+Node_radius*(2.1+((double)rand()/(double)RAND_MAX)*0.2)*sin(theta)*sin(phi);
//        double temp_z=Node_Position[i-1][2]+Node_radius*(2.1+((double)rand()/(double)RAND_MAX)*0.2)*cos(theta);
        double temp_x=Node_Position[i-1][0]+Node_radius*2.1*sin(theta)*cos(phi);
        double temp_y=Node_Position[i-1][1]+Node_radius*2.1*sin(theta)*sin(phi);
        double temp_z=Node_Position[i-1][2]+Node_radius*2.1*cos(theta);
        
        bool accept_random_step=true;
        for (int k=0; k<i; k++) {
            double temp_dist=sqrt( (temp_x-Node_Position[k][0])*(temp_x-Node_Position[k][0]) + (temp_y-Node_Position[k][1])*(temp_y-Node_Position[k][1]) + (temp_z-Node_Position[k][2])*(temp_z-Node_Position[k][2])  );
            if (temp_dist<Node_radius*2.1){
                i--;
                accept_random_step=false;
                break;
            }
        } // for (int k=0; k<i; k++)
        
        if (accept_random_step) {
            Node_Position[i][0] = temp_x;
            Node_Position[i][1] = temp_y;
            Node_Position[i][2] = temp_z;
//            cout<<"COM "<<i<<": \n"<<temp_x<<"\t"<<temp_y<<"\t"<<temp_z<<endl;
            
            Node_Velocity[i][0]=((double)rand()/(double)RAND_MAX)*2-1;
            Node_Velocity[i][1]=((double)rand()/(double)RAND_MAX)*2-1;
            Node_Velocity[i][2]=((double)rand()/(double)RAND_MAX)*2-1;
            
            position_COM[0] += temp_x;
            position_COM[1] += temp_y;
            position_COM[2] += temp_z;

            velocity_COM[0]+=Node_Velocity[i][0];
            velocity_COM[1]+=Node_Velocity[i][1];
            velocity_COM[2]+=Node_Velocity[i][2];
        }
        
    } // for (int i=1; i<Num_of_Nodes; i++)
    
    
    
    position_COM[0] /= Num_of_Nodes;
    position_COM[1] /= Num_of_Nodes;
    position_COM[2] /= Num_of_Nodes;
    
    velocity_COM[0]/=Num_of_Nodes;
    velocity_COM[1]/=Num_of_Nodes;
    velocity_COM[2]/=Num_of_Nodes;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Position[i][0] -= position_COM[0];
        Node_Position[i][1] -= position_COM[1];
        Node_Position[i][2] -= position_COM[2];
        
        Node_Velocity[i][0] -= velocity_COM[0];
        Node_Velocity[i][1] -= velocity_COM[1];
        Node_Velocity[i][2] -= velocity_COM[2];
    }
    
    Membrane_neighbour_node.resize(Num_of_Nodes);
}
