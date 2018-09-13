//
//  Membrane_node_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Node_Bonds_identifier(void){
    
    vector<int> push;
    push.resize(2);
    
//    int temp_Membrane_num_of_Node_Pairs=0;
    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Num_of_Triangles;i++)
    {
        temp_Membrane_triangle_Node_A= Triangle_list[i][0];
        temp_Membrane_triangle_Node_B= Triangle_list[i][1];
        temp_Membrane_triangle_Node_C= Triangle_list[i][2];
        
        for(int j=0;j<Node_Bond_list.size();j++)
        {
            if(  ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_A &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_B )  || ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_B &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_B &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_C &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_A &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Node_Bond_list[j][0]==temp_Membrane_triangle_Node_C &  Node_Bond_list[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            push[0]=temp_Membrane_triangle_Node_A;
            push[1]=temp_Membrane_triangle_Node_B;
            Node_Bond_list.push_back(push);
        }
        
        if(repeatednumber2==0)
        {
            push[0]=temp_Membrane_triangle_Node_B;
            push[1]=temp_Membrane_triangle_Node_C;
            Node_Bond_list.push_back(push);
            
        }
        
        if(repeatednumber3==0)
        {
            push[0]=temp_Membrane_triangle_Node_A;
            push[1]=temp_Membrane_triangle_Node_C;
            Node_Bond_list.push_back(push);
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    Num_of_Node_Pairs=Node_Bond_list.size();
    cout<<"Membrane # of node pairs: "<<Num_of_Node_Pairs<<endl;
}

