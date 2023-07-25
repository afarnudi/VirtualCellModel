#include "Membrane.h"
//i guess there is no need to define this function. we know that for each triangle, there are 3 other triangles which share an edge with the first triangle, so we have Membrane_num_of_Triangle_Pairs =3*(membrane_trianle_list.size())/2. in fact the results of this function for different meshes are the same as this simple calculations. but i have a question about the mesh. we know that it is impossible to cover a spherical shell with triangles. how does gmesh do this? and in the cases other than spherical shell (like RBC) does this calculation work or not? 

void Membrane::Triangle_pair_counter ()
{
    
        //In this function we count the total number of triangles that have a common edge (we count them twice, hence report half the number at the end).
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C;
    Num_of_Triangle_Pairs=0;  // This counts the number of triangle pairs that have an edge in common.
    for(int i=0 ;i<Num_of_Triangles-1;i++)  // who are neighbors??
    {
        temp_triangle_node_A=Triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        temp_triangle_node_B=Triangle_list[i][1];
        temp_triangle_node_C=Triangle_list[i][2];
        
        for(int j=i+1;j<Num_of_Triangles;j++)
        {
            if      ( (Triangle_list[j][0]==temp_triangle_node_A)  &&  (Triangle_list[j][1]==temp_triangle_node_B)  && (Triangle_list[j][2]!=temp_triangle_node_C) ){
                Num_of_Triangle_Pairs++;
            }
            else if     ( (Triangle_list[j][0]==temp_triangle_node_B)  &&  (Triangle_list[j][1]==temp_triangle_node_A)  && (Triangle_list[j][2]!=temp_triangle_node_C) ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]!=temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]!=temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            else if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]!=temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            else if     ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]!=temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_B  &&  Triangle_list[j][1]!=temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_C ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]!=temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_C ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            
            else if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]!=temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            else if     ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]!=temp_triangle_node_B ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_C  &&  Triangle_list[j][1]!=temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]==temp_triangle_node_A  &&  Triangle_list[j][1]!=temp_triangle_node_B  && Triangle_list[j][2]==temp_triangle_node_C ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_C  && Triangle_list[j][2]==temp_triangle_node_A ){
                Num_of_Triangle_Pairs++;
            }
            else if      ( Triangle_list[j][0]!=temp_triangle_node_B  &&  Triangle_list[j][1]==temp_triangle_node_A  && Triangle_list[j][2]==temp_triangle_node_C ){
                Num_of_Triangle_Pairs++;
            }
        }
    }
//    Num_of_Triangle_Pairs=Num_of_Triangle_Pairs/2;
}
