//
//  Membrane_node_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"

using namespace std;

void add_bond_to_list (vector<vector<int> > tempbonds, vector<vector<int> > &Node_Bond_list){
    for (int i=0; i<tempbonds.size(); i++) {
        bool add=true;
        for (int j=0; j<Node_Bond_list.size(); j++) {
            if (tempbonds[i][0]==Node_Bond_list[j][0] && tempbonds[i][1]==Node_Bond_list[j][1]) {
                add=false;
                break;
            }
            if (tempbonds[i][0]==Node_Bond_list[j][1] && tempbonds[i][1]==Node_Bond_list[j][0]) {
                add=false;
                break;
            }
            
        }
        if (add) {
            Node_Bond_list.push_back(tempbonds[i]);
        }
    }
}

void Actin::Node_Bond_identifier_3(void){
    
    Node_Bond_list.clear();
    for(int i=0;i<Pyramid_Nodes.size();i++)
    {
        vector<vector<int> > tempbonds;
        vector<int> push;
        push.resize(2);
        push[0]=Pyramid_Nodes[i][0];
        push[1]=Pyramid_Nodes[i][1];
        tempbonds.push_back(push);
        push[0]=Pyramid_Nodes[i][0];
        push[1]=Pyramid_Nodes[i][2];
        tempbonds.push_back(push);
        push[0]=Pyramid_Nodes[i][0];
        push[1]=Pyramid_Nodes[i][3];
        tempbonds.push_back(push);
        push[0]=Pyramid_Nodes[i][1];
        push[1]=Pyramid_Nodes[i][2];
        tempbonds.push_back(push);
        push[0]=Pyramid_Nodes[i][1];
        push[1]=Pyramid_Nodes[i][3];
        tempbonds.push_back(push);
        push[0]=Pyramid_Nodes[i][2];
        push[1]=Pyramid_Nodes[i][3];
        tempbonds.push_back(push);
        
        add_bond_to_list(tempbonds,Node_Bond_list);
        
    }
    
    Num_of_Node_Pairs=int(Node_Bond_list.size());
    cout<<"# of node pairs: "<<Num_of_Node_Pairs<<endl;
}


void Actin::Node_Bond_identifier(void){
    
    vector<int> push;
    push.resize(2);
    
    //The first Pyramid:
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][1];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][2];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][1];
    push[1]=Pyramid_Nodes[0][2];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][1];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][2];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);
    
    int pyramid_Node_1, pyramid_Node_2, pyramid_Node_3, pyramid_Node_4;
    
    bool repeated_pair_1=false;
    bool repeated_pair_2=false;
    bool repeated_pair_3=false;
    bool repeated_pair_4=false;
    bool repeated_pair_5=false;
    bool repeated_pair_6=false;
    
    
    
    for(int i=1;i<Pyramid_Nodes.size();i++)
    {
        pyramid_Node_1=Pyramid_Nodes[i][0];
        pyramid_Node_2=Pyramid_Nodes[i][1];
        pyramid_Node_3=Pyramid_Nodes[i][2];
        pyramid_Node_4=Pyramid_Nodes[i][3];
        
        for(int j=0;j<Node_Bond_list.size();j++)
        {
            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_2 )  || (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_1=true;
            }
            
            if(  (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_3 )  || (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_2 )    )
            {
                repeated_pair_2=true;
            }
            
            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_3 )  || (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_3=true;
            }
            
            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_4=true;
            }
            
            if(  (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_2 )    )
            {
                repeated_pair_5=true;
            }
            
            if(  (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_3 )    )
            {
                repeated_pair_6=true;
            }
        }
        
        if(!repeated_pair_1)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_2;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_2)
        {
            push[0]=pyramid_Node_2;
            push[1]=pyramid_Node_3;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_3)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_3;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_4)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_5)
        {
            push[0]=pyramid_Node_2;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_6)
        {
            push[0]=pyramid_Node_3;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }
        
        repeated_pair_1=false;
        repeated_pair_2=false;
        repeated_pair_3=false;
        repeated_pair_4=false;
        repeated_pair_5=false;
        repeated_pair_6=false;
        
    }
    
    Num_of_Node_Pairs=int(Node_Bond_list.size());
//    cout<<"# of node pairs: "<<Num_of_Node_Pairs<<endl;
}

// for test
void Actin::Node_Bond_identifier_2(void){
    for(int i=0; i<filaments.size(); i++)
    {
        Node_Bond_list.push_back(filaments[i]);
        //std::cout<<filaments[i][0]<<"and"<<filaments[i][1]<<'\n';
    }
    Num_of_Node_Pairs=int(Node_Bond_list.size());
    
    for(int i=0; i<abps.size(); i++)
    {
        abp_Bond_list.push_back(abps[i]);
        //std::cout<<filaments[i][0]<<"and"<<filaments[i][1]<<'\n';
    }
    Num_of_abp_Pairs=int(abp_Bond_list.size());
    
    
    for(int i=0; i<MTs.size(); i++)
    {
        MT_Bond_list.push_back(MTs[i]);
        //std::cout<<filaments[i][0]<<"and"<<filaments[i][1]<<'\n';
    }
    Num_of_MT_Pairs=int(MT_Bond_list.size());
}

