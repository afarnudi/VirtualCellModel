//
//  Membrane_triangle_pair_identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include <stdio.h>
#include "Membrane.h"

void Membrane::Triangle_pair_identifier(void){
    Triangle_Pair_Nodes.clear();
    Triangle_pair_list.clear();
    Triangle_Pair_Nodes.resize(Num_of_Triangle_Pairs);
    
    vector<int> temp_triangle_pair;
    temp_triangle_pair.resize(2);
    
    int temp_int_Triangle_Pair_index=0;
    
//    int triangle_pairs=0;
//    int temp[4][2*Membrane_num_of_Triangle_Pairs];
//    int temp2[4][2*Membrane_num_of_Triangle_Pairs];
    
    
    for(int i=0 ;i<Num_of_Triangles;i++)
    {
        int neighbour=-1;
        int temp_triangle_node_A=Triangle_list[i][0];
        int temp_triangle_node_B=Triangle_list[i][1];
        int temp_triangle_node_C=Triangle_list[i][2];
        int neighbour_indicator=0;
        //        Indicates the existence of Node neighbour for a node pair (other than the other node of triangle):
        //        neighbour_indicator=0 No Node Pairs; neighbour_indicator=1, for temp_triangle_node_A-temp_triangle_node_B; neighbour_indicator=2, for temp_triangle_node_B-temp_triangle_node_C; And neighbour_indicator=3 for temp_triangle_node_C-temp_triangle_node_A.
        for(int j=i+1;j<Num_of_Triangles;j++)
        {
            //************************** finding neighbours **************************
            // neibours of temp_triangle_node_A-temp_triangle_node_B:
            if ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]!=temp_triangle_node_C ){
                neighbour=Triangle_list[j][2];
                neighbour_indicator=1;
            }
            if     ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]!=temp_triangle_node_C )
            {
                neighbour=Triangle_list[j][2];
                neighbour_indicator=1;
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]!=temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]!=temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=1;
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=1;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Triangle_list[j][2];
                neighbour_indicator=2;
            }
            if     ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Triangle_list[j][2];
                neighbour_indicator=2;
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]!=temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Triangle_list[j][1];
                neighbour_indicator=2;
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]!=temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_B )
            {
                neighbour=Triangle_list[j][1];
                neighbour_indicator=2;
                
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=2;
                
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=2;
            }
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if     ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]!=temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]!=temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=3;
            }
            if      ( Triangle_list[j][0]!=temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Triangle_list[j][0];
                neighbour_indicator=3;
            }
            
            if(neighbour_indicator!=0)  //  to speed up  the programme we first check if we have found a neighbour or not
            {
                // note that temp_triangle_node_A-temp_triangle_node_B-temp_triangle_node_C-neighbour  are 4 point of two triangle wich will interact
                Triangle_Pair_Nodes[temp_int_Triangle_Pair_index].resize(4);
                if(neighbour_indicator==1)
                {
                    
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_C;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_A;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_B;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                    temp_triangle_pair[0]=i;
                    temp_triangle_pair[1]=j;
                    Triangle_pair_list.push_back(temp_triangle_pair);
                    
                    
                } else if(neighbour_indicator==2)
                {
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_A;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_C;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_B;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                    temp_triangle_pair[0]=i;
                    temp_triangle_pair[1]=j;
                    Triangle_pair_list.push_back(temp_triangle_pair);
                    
                } else if (neighbour_indicator==3)
                {
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][0]=temp_triangle_node_B;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][1]=temp_triangle_node_C;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][2]=temp_triangle_node_A;
                    Triangle_Pair_Nodes[temp_int_Triangle_Pair_index][3]=neighbour;
                    temp_triangle_pair[0]=i;
                    temp_triangle_pair[1]=j;
                    Triangle_pair_list.push_back(temp_triangle_pair);
                }
                temp_int_Triangle_Pair_index++;
            }
            neighbour_indicator=0;
        }
    }
    if (Triangle_pair_list.size()!=Num_of_Triangle_Pairs) {
        cout<<"Triangle_pair_list.size()!=Num_of_Triangle_Pairs\nTriangle_pair_list.size() "<<Triangle_pair_list.size()<<"\nNum_of_Triangle_Pairs "<<Num_of_Triangle_Pairs<<endl;
    }
    
    
}
