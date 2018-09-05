//
//  ECM_Node_Pair_Identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.hpp"

void ECM::ECM_Node_Pair_Identifier(void){
    ECM_num_of_Node_Pairs=0;
    vector<vector<int> > ECM_Node_Pair_list;
    vector<int> Node_Pairs;
    Node_Pairs.resize(2);
    
//    int ECM_pyramid_Node_1, ECM_pyramid_Node_2, ECM_pyramid_Node_3, ECM_pyramid_Node_4;
//    double temp_double;
//    int Num_of_objects;
    
    int temp_ECM_triangle_Node_A, temp_ECM_triangle_Node_B, temp_ECM_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<ECM_num_of_Triangles;i++)
    {
        temp_ECM_triangle_Node_A= ECM_triangle_list[i][0];
        temp_ECM_triangle_Node_B= ECM_triangle_list[i][1];
        temp_ECM_triangle_Node_C= ECM_triangle_list[i][2];
        
        for(int j=0;j<ECM_Node_Pair_list.size();j++)
        {
            if(  ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_A &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_B )  || ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_B &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_B &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_C )  || ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_C &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_A &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_C )  || ( ECM_Node_Pair_list[j][0]==temp_ECM_triangle_Node_C &&  ECM_Node_Pair_list[j][1]==temp_ECM_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            Node_Pairs[0]=temp_ECM_triangle_Node_A;
            Node_Pairs[1]=temp_ECM_triangle_Node_B;
            
            ECM_Node_Pair_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber2==0)
        {
            Node_Pairs[0]=temp_ECM_triangle_Node_B;
            Node_Pairs[1]=temp_ECM_triangle_Node_C;
            
            ECM_Node_Pair_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber3==0)
        {
            Node_Pairs[0]=temp_ECM_triangle_Node_A;
            Node_Pairs[1]=temp_ECM_triangle_Node_C;
            
            ECM_Node_Pair_list.push_back(Node_Pairs);
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    ECM_num_of_Node_Pairs=ECM_Node_Pair_list.size();
    cout<<"ECM # of node pairs: "<<ECM_Node_Pair_list.size()<<endl;
}

