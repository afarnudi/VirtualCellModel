//
//  ECM_Node_Pair_Identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"

void ECM::Node_Bond_identifier(void){
    
    Num_of_Node_Pairs=0;
    
    vector<int> Node_Pairs;
    Node_Pairs.resize(2);
    
//    int ECM_pyramid_Node_1, ECM_pyramid_Node_2, ECM_pyramid_Node_3, ECM_pyramid_Node_4;
//    double temp_double;
//    int Num_of_objects;
    
    int triangle_Node_A, triangle_Node_B, triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Num_of_Triangles;i++)
    {
        triangle_Node_A= Triangle_List[i][0];
        triangle_Node_B= Triangle_List[i][1];
        triangle_Node_C= Triangle_List[i][2];
        
        for(int j=0;j<Node_Pair_list.size();j++)
        {
            if(  ( Node_Pair_list[j][0]==triangle_Node_A &&  Node_Pair_list[j][1]==triangle_Node_B )  ||
                 ( Node_Pair_list[j][0]==triangle_Node_B &&  Node_Pair_list[j][1]==triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( Node_Pair_list[j][0]==triangle_Node_B &&  Node_Pair_list[j][1]==triangle_Node_C )  ||
                 ( Node_Pair_list[j][0]==triangle_Node_C &&  Node_Pair_list[j][1]==triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( Node_Pair_list[j][0]==triangle_Node_A &&  Node_Pair_list[j][1]==triangle_Node_C )  ||
                 ( Node_Pair_list[j][0]==triangle_Node_C &&  Node_Pair_list[j][1]==triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            Node_Pairs[0]=triangle_Node_A;
            Node_Pairs[1]=triangle_Node_B;
            
            Node_Pair_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber2==0)
        {
            Node_Pairs[0]=triangle_Node_B;
            Node_Pairs[1]=triangle_Node_C;
            
            Node_Pair_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber3==0)
        {
            Node_Pairs[0]=triangle_Node_A;
            Node_Pairs[1]=triangle_Node_C;
            
            Node_Pair_list.push_back(Node_Pairs);
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    Num_of_Node_Pairs=int(Node_Pair_list.size());
    cout<<"ECM # of node pairs: "<<Node_Pair_list.size()<<endl;
}

