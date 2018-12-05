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
    cout<<"Initialising the Chromatin Class..."<<endl;
    build_random_chain();
    cout<<"\n\nChromatinb class initiated.\n";

}

void Chromatin::build_random_chain(void){
    
    Node_Position.resize(Num_of_Nodes);
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
    }
    Node_Position[0][0]=0;
    Node_Position[0][1]=0;
    Node_Position[0][2]=0;
    
    const double PI  =3.141592653589793238463;
    
    for (int i=1; i<Num_of_Nodes; i++) {
        double theta=((double)rand()/(double)RAND_MAX)*PI;
        double phi=((double)rand()/(double)RAND_MAX)*2*PI;
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
            Node_Position[i][0]=temp_x;
            Node_Position[i][1]=temp_y;
            Node_Position[i][2]=temp_z;
        }
        
    } // for (int i=1; i<Num_of_Nodes; i++)
    
}
