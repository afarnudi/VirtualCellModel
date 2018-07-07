#include "Membrane.h"
//i guess there is no need to define this function. we know that for each triangle, there are 3 other triangles which share an edge with the first triangle, so we have Membrane_num_of_Triangle_Pairs =3*(membrane_trianle_list.size())/2. in fact the results of this function for different meshes are the same as this simple calculations. but i have a question about the mesh. we know that it is impossible to cover a spherical shell with triangles. how does gmesh do this? and in the cases other than spherical shell (like RBC) does this calculation work or not? 

int Membrane::Membrane_triangle_pair_counter ()
{
	
	    //In this function we count the total number of triangles that have a common edge (we count them twice, hence report half the number at the end).
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C;
    int triangle_pairs=0;  // This counts the number of triangle pairs that have an edge in common.
    for(int i=0 ;i<Membrane_triangle_list.size();i++)  // who are neighbors??
    {
        temp_triangle_node_A=Membrane_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        temp_triangle_node_B=Membrane_triangle_list[i][1];
        temp_triangle_node_C=Membrane_triangle_list[i][2];
        
        for(int j=0;j<Membrane_triangle_list.size();j++)
        {
            if      ( (Membrane_triangle_list[j][0]==temp_triangle_node_A)  &  (Membrane_triangle_list[j][1]==temp_triangle_node_B)  & (Membrane_triangle_list[j][2]!=temp_triangle_node_C) ){
                triangle_pairs++;
            }
            if     ( (Membrane_triangle_list[j][0]==temp_triangle_node_B)  &  (Membrane_triangle_list[j][1]==temp_triangle_node_A)  & (Membrane_triangle_list[j][2]!=temp_triangle_node_C) ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                triangle_pairs++;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                triangle_pairs++;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
        }
    }
    return (triangle_pairs/2);
}