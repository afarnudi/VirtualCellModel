//
//  Initialise.cpp
//  Membrae
//
//  Created by Ali Farnudi on 14/10/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"
#include <algorithm>
//void build_random_chain(void);

using std::cout;
using std::endl;



void Chromatin::build_random_chain(void){
    
    
    
    Node_Velocity.resize(Num_of_Nodes);
    Node_Force.resize(Num_of_Nodes);
    
    
    double velocity_COM[3]={0};
    
    for(int i=0;i<Num_of_Nodes;i++)
    {
        Node_Velocity[i].resize(3,0);
        Node_Force[i].resize(3,0);
    }
//    const double PI  =3.141592653589793238463;
    int attempt=0;
    cout<<endl;
    try{
      random_walk_gen(velocity_COM);
    }
    catch (int e){
      
        attempt++;
//        cout<<"attempt "<<attempt<<endl;
        printf("attempt # %d \r",attempt);
        cout<< std::flush;
        cout<<"attempt "<<attempt<<endl;
        if (attempt>2000) {
            cout<<endl;
            cout<<"I have tried and I can't generate a self avoiding chain with these parameters. \n";
            exit(0);
        }
        random_walk_gen(velocity_COM);
    }
    cout<<endl;
    velocity_COM[0]/=Num_of_Nodes;
    velocity_COM[1]/=Num_of_Nodes;
    velocity_COM[2]/=Num_of_Nodes;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Velocity[i][j] -= velocity_COM[j];
        }
    }
    
    cout<<"Generated a self avoiding random walk chain."<<endl;
}

void Chromatin::random_walk_gen(double velocity_COM[3]){
    
    velocity_COM[0]=0;
    velocity_COM[1]=0;
    velocity_COM[2]=0;
    Node_Position.clear();
    Node_Position.resize(Num_of_Nodes);
    for(auto &pos:Node_Position){
        pos.resize(3,0);
    }
    int attempt_counter=0;
    srand(time(NULL));
    
    double bond_length = stod(BondNominalLength_stat);
    
    for (int i=1; i<Num_of_Nodes; i++) {
        double theta=((double)rand()/(double)RAND_MAX)*M_PI;
        double phi=((double)rand()/(double)RAND_MAX)*2*M_PI;
        double temp_x=Node_Position[i-1][0]+bond_length*sin(theta)*cos(phi);
        double temp_y=Node_Position[i-1][1]+bond_length*sin(theta)*sin(phi);
        double temp_z=Node_Position[i-1][2]+bond_length*cos(theta);
        
        bool accept_random_step=true;
        for (int k=0; k<i-1; k++) {
            double temp_dist=sqrt( (temp_x-Node_Position[k][0])*(temp_x-Node_Position[k][0]) + (temp_y-Node_Position[k][1])*(temp_y-Node_Position[k][1]) + (temp_z-Node_Position[k][2])*(temp_z-Node_Position[k][2])  );
            if (temp_dist<Node_radius*2.1){
                i--;
                accept_random_step=false;
                attempt_counter++;
                if (attempt_counter>2000) {
                    string errorMessage = TWARN;
                    errorMessage+="Chroamtin coordinate generator: Self Avoiding Random Walk: "+std::to_string(attempt_counter)+"attempts failed in generating a SAW chain. Please try a shorter chain or the seed.";
                    errorMessage+= TRESET;
                    throw std::runtime_error(errorMessage);
                }
                break;
            }
        } // for (int k=0; k<i; k++)
        
        if (accept_random_step) {
            attempt_counter=0;
            
            Node_Position[i][0] = temp_x;
            Node_Position[i][1] = temp_y;
            Node_Position[i][2] = temp_z;
            
//            Node_Velocity[i][0]=((double)rand()/(double)RAND_MAX)*2-1;
//            Node_Velocity[i][1]=((double)rand()/(double)RAND_MAX)*2-1;
//            Node_Velocity[i][2]=((double)rand()/(double)RAND_MAX)*2-1;
//            velocity_COM[0] += Node_Velocity[i][0];
//            velocity_COM[1] += Node_Velocity[i][1];
//            velocity_COM[2] += Node_Velocity[i][2];
        }
    } // for (int i=1; i<Num_of_Nodes; i++)
}
