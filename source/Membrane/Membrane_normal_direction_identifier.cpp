#include "Membrane.h"
#include "General_functions.hpp"

void Membrane::Normal_direction_Identifier()
{
    double AC[3], AB[3], ABxAC[3], reference[3];
    int Point_A, Point_B, Point_C;
    
    
    for(  int i=0;i<Triangle_list.size();i++)
    {
        Point_A=Triangle_list[i][0];
        Point_B=Triangle_list[i][1];
        Point_C=Triangle_list[i][2];
        
        for (int index =0; index<3; index++) {
            AB[index]=Node_Position[Point_B][index]-Node_Position[Point_A][index];
            AC[index]=Node_Position[Point_C][index]-Node_Position[Point_A][index];
            if (UseXYZinMembrane) {
                reference[index]=XYZinMembrane[index]-Node_Position[Point_A][index];
            } else {
                reference[index]=COM_position[index]-Node_Position[Point_A][index];
            }
        }
        
        crossvector(ABxAC, AB, AC);
        //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Triangle_list[][0], B=Triangle_list[][1], C=Triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
        
        if(innerproduct(ABxAC, reference)<0 )
        {
            Triangle_list[i][1]=Point_C;
            Triangle_list[i][2]=Point_B;
            //cout<<"min"<<endl;
        }
    } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
}// END OF: Normal_direction_Identifier function
//end of normal_direction_identifire***
