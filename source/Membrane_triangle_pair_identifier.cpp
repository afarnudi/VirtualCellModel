//
//  Membrane_triangle_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include <stdio.h>
#include "Membrane.h"

void Membrane::membrane_triangle_pair_identifier(void){
    
    Membrane_Triangle_Pair_Nodes.resize(Num_of_Triangle_Pairs);
    int temp_int_Triangle_Pair_index=0;
    
//    int triangle_pairs=0;
//    int temp[4][2*Membrane_num_of_Triangle_Pairs];
//    int temp2[4][2*Membrane_num_of_Triangle_Pairs];
    
    
    for(int i=0 ;i<Membrane_num_of_Triangles;i++)
    {
        int neighbour=-1;
        int temp_triangle_node_A=Membrane_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        int temp_triangle_node_B=Membrane_triangle_list[i][1];
        int temp_triangle_node_C=Membrane_triangle_list[i][2];
        int neighbour_indicator=0;
        //        Indicates the existence of Node neighbour for a node pair (other than the membrane of the triangle):
        //        neighbour_indicator=0 No Node Pairs; neighbour_indicator=1, for temp_triangle_node_A-temp_triangle_node_B; neighbour_indicator=2, for temp_triangle_node_B-temp_triangle_node_C; And neighbour_indicator=3 for temp_triangle_node_C-temp_triangle_node_A.
        for(int j=i;j<Membrane_num_of_Triangles;j++)
        {
            //************************** finding neighbours **************************
            // neibours of temp_triangle_node_A-temp_triangle_node_B:
            if ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=1;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_C )
            {
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=1;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=2;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B )
            {
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=2;
            }
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=3;
            }
            
            if(neighbour_indicator!=0)  //  to speed up  the programme we first check if we have found a neighbour or not
            {
                // note that temp_triangle_node_A-temp_triangle_node_B-temp_triangle_node_C-neighbour  are 4 point of two triangle wich will interact
                Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index].resize(4);
                if(neighbour_indicator==1)
                {
                    
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_C;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_A;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_B;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                    
                    
                } else if(neighbour_indicator==2)
                {
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_A;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_C;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_B;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                    
                } else if (neighbour_indicator==3)
                {
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_B;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_C;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_A;
                    Membrane_Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                }
                temp_int_Triangle_Pair_index++;
            }
            neighbour_indicator=0;
        }
    }
    
}
