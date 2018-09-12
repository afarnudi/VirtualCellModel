//
//  Membrane_node_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::membrane_node_pair_identifier(void){
    
    vector<int> push;
    push.resize(2);
    
//    int temp_Membrane_num_of_Node_Pairs=0;
    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Membrane_num_of_Triangles;i++)
    {
        temp_Membrane_triangle_Node_A= Membrane_triangle_list[i][0];
        temp_Membrane_triangle_Node_B= Membrane_triangle_list[i][1];
        temp_Membrane_triangle_Node_C= Membrane_triangle_list[i][2];
        
        for(int j=0;j<Membrane_Edges.size();j++)
        {
            if(  ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_B )  || ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Edges[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Edges[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            push[0]=temp_Membrane_triangle_Node_A;
            push[1]=temp_Membrane_triangle_Node_B;
            Membrane_Edges.push_back(push);
        }
        
        if(repeatednumber2==0)
        {
            push[0]=temp_Membrane_triangle_Node_B;
            push[1]=temp_Membrane_triangle_Node_C;
            Membrane_Edges.push_back(push);
            
        }
        
        if(repeatednumber3==0)
        {
            push[0]=temp_Membrane_triangle_Node_A;
            push[1]=temp_Membrane_triangle_Node_C;
            Membrane_Edges.push_back(push);
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    Num_of_Node_Pairs=Membrane_Edges.size();
    cout<<"Membrane # of node pairs: "<<Num_of_Node_Pairs<<endl;
}

