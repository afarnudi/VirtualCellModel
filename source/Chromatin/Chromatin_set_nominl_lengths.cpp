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

using namespace std;

void Chromatin::set_bond_nominal_lengths(void){
    BondNominalLengths.clear();
    BondNominalLengths.resize(Num_of_Nodes-1,0);
    
    Min_node_pair_length=100000000000000000000.;
    Max_node_pair_length=0;
    Average_node_pair_length=0;
    
    
    for (int i=1; i<Num_of_Nodes; i++) {
        double bondvec[3]={Node_Position[i][0]-Node_Position[i-1][0],
                           Node_Position[i][1]-Node_Position[i-1][1],
                           Node_Position[i][2]-Node_Position[i-1][2]};
        BondNominalLengths[i-1]=vector_length(bondvec);
        
        Average_node_pair_length += BondNominalLengths[i-1];
        
        if (BondNominalLengths[i-1] < Min_node_pair_length) {
            Min_node_pair_length = BondNominalLengths[i-1];
        }
        if (BondNominalLengths[i-1] > Max_node_pair_length) {
            Max_node_pair_length = BondNominalLengths[i-1];
        }
        
    }
    Average_node_pair_length/=Num_of_Nodes-1;
    if (BondNominalLength_stat=="Au") {
        if (!generalParameters.Testmode) {
            cout<<"Node pair (bond) distances:\n";
            cout<<"\tMax "<<Max_node_pair_length<<"\tMin "<<Min_node_pair_length<<"\tAverage "<<Average_node_pair_length<<endl;
        }
    } else if (BondNominalLength_stat=="Av"){
        BondNominalLengths.clear();
        BondNominalLengths.resize(Num_of_Nodes-1,Average_node_pair_length);
        if (!generalParameters.Testmode) {
            cout<<"Average node pair (bond) distance "<<TFILE<<Average_node_pair_length<<TRESET<<" will be set for the spring nominal length."<<endl;
        }
    } else {
        BondNominalLengths.clear();
        BondNominalLengths.resize(Num_of_Nodes,stod(BondNominalLength_stat));
        if (!generalParameters.Testmode) {
            cout<<"User defined nominal length ("<<TFILE<<BondNominalLength_stat<<TRESET<<") will be used."<<endl;
        }
    }
    
}

void Chromatin::set_node_radius(void){
    
    if (Node_radius_stat=="Av") {
        Node_radius=0.5*Average_node_pair_length;
        cout<<"Node radius set to half of bond average distances."<<endl;
    } else {
        Node_radius=stod(Node_radius_stat);
    }
}

void Chromatin::SpontaneousBendingAngle(void){
    if (Num_of_Nodes<3) {
        string errorMessage = TWARN;
        errorMessage +="Chromatin initliser: The chromatin must have at least 3 nodes in order to set up an angle potential. You can disable the angle potential or use a longer chromatin.";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    SpontaneousBendingAngleInRad = SpontaneousBendingAngleInDeg*M_PI/180.;
}
