//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

using std::cout;
using std::endl;

void Actin::check(void){
    Min_node_pair_length=1000;
    Max_node_pair_length=0;
    Average_node_pair_length=0;
    for (int i=0; i<Num_of_Node_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])*(Node_Position[Node_Bond_list[i][0]][0]-Node_Position[Node_Bond_list[i][1]][0])+(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])*(Node_Position[Node_Bond_list[i][0]][1]-Node_Position[Node_Bond_list[i][1]][1])+(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2])*(Node_Position[Node_Bond_list[i][0]][2]-Node_Position[Node_Bond_list[i][1]][2]));
        Average_node_pair_length+=dist;
        
        if (Min_node_pair_length>dist) {
            Min_node_pair_length=dist;
        }
        if (dist>Max_node_pair_length) {
            Max_node_pair_length=dist;
        }
    }
    Average_node_pair_length/=Num_of_Node_Pairs;
    
    if (!generalParameters.Testmode) {
        cout<<"Node pair (bond) distances:\n";
        cout<<"\tMax "<<Max_node_pair_length<<"\tMin "<<Min_node_pair_length<<"\tAverage "<<Average_node_pair_length<<endl;
    }
}


void Actin::check_2(void){
    Min_abp_pair_length=1000;
    Max_abp_pair_length=0;
    Average_abp_pair_length=0;
    for (int i=0; i<Num_of_abp_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[abp_Bond_list[i][0]][0]-Node_Position[abp_Bond_list[i][1]][0])*(Node_Position[abp_Bond_list[i][0]][0]-Node_Position[abp_Bond_list[i][1]][0])+(Node_Position[abp_Bond_list[i][0]][1]-Node_Position[abp_Bond_list[i][1]][1])*(Node_Position[abp_Bond_list[i][0]][1]-Node_Position[abp_Bond_list[i][1]][1])+(Node_Position[abp_Bond_list[i][0]][2]-Node_Position[abp_Bond_list[i][1]][2])*(Node_Position[abp_Bond_list[i][0]][2]-Node_Position[abp_Bond_list[i][1]][2]));
        Average_abp_pair_length+=dist;
        
        if (Min_abp_pair_length>dist) {
            Min_abp_pair_length=dist;
        }
        if (dist>Max_abp_pair_length) {
            Max_abp_pair_length=dist;
        }
    }
    Average_abp_pair_length/=Num_of_abp_Pairs;
    
    cout<<"Max abp distance="<<Max_abp_pair_length<<"\tmin abp distance="<<Min_abp_pair_length<<"\tAverage abp distance="<<Average_abp_pair_length<<endl;

}








void Actin::check_3(void){
    Min_MT_pair_length=1000;
    Max_MT_pair_length=0;
    Average_MT_pair_length=0;
    for (int i=0; i<Num_of_MT_Pairs; i++) {
        double dist=0;
        dist=sqrt((Node_Position[MT_Bond_list[i][0]][0]-Node_Position[MT_Bond_list[i][1]][0])*(Node_Position[MT_Bond_list[i][0]][0]-Node_Position[MT_Bond_list[i][1]][0])+(Node_Position[MT_Bond_list[i][0]][1]-Node_Position[MT_Bond_list[i][1]][1])*(Node_Position[MT_Bond_list[i][0]][1]-Node_Position[MT_Bond_list[i][1]][1])+(Node_Position[MT_Bond_list[i][0]][2]-Node_Position[MT_Bond_list[i][1]][2])*(Node_Position[MT_Bond_list[i][0]][2]-Node_Position[MT_Bond_list[i][1]][2]));
        Average_MT_pair_length+=dist;
        
        if (Min_MT_pair_length>dist) {
            Min_MT_pair_length=dist;
        }
        if (dist>Max_MT_pair_length) {
            Max_MT_pair_length=dist;
        }
    }
    Average_MT_pair_length/=Num_of_MT_Pairs;
    
    cout<<"Max MT distance="<<Max_MT_pair_length<<"\tmin MT distance="<<Min_MT_pair_length<<"\tAverage MT distance="<<Average_MT_pair_length<<endl;

}


