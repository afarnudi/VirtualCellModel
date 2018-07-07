#include "Membrane.h"
#include "General_functions.hpp"

void Membrane::Membrane_Normal_direction_Identifier()
{
    double AC[3], AB[3], ABxAC[3], xyz[3];
    int Point_A, Point_B, Point_C;
    cout<<Membrane_triangle_list.size()<<endl;
    for(  int i=0;i<Membrane_triangle_list.size();i++)
    {
        Point_A=Membrane_triangle_list[i][0];
        Point_B=Membrane_triangle_list[i][1];
        Point_C=Membrane_triangle_list[i][2];
        
        AB[0]=Membrane_Node_Position[Point_B][0]-Membrane_Node_Position[Point_A][0];
        AB[1]=Membrane_Node_Position[Point_B][1]-Membrane_Node_Position[Point_A][1];
        AB[2]=Membrane_Node_Position[Point_B][2]-Membrane_Node_Position[Point_A][2];
        
        AC[0]=Membrane_Node_Position[Point_C][0]-Membrane_Node_Position[Point_A][0];
        AC[1]=Membrane_Node_Position[Point_C][1]-Membrane_Node_Position[Point_A][1];
        AC[2]=Membrane_Node_Position[Point_C][2]-Membrane_Node_Position[Point_A][2];
        
        xyz[0]=Membrane_Node_Position[Point_A][0];
        xyz[1]=Membrane_Node_Position[Point_A][1];
        xyz[2]=Membrane_Node_Position[Point_A][2];
        
        crossvector(ABxAC, AB, AC);
        //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
        
        if(innerproduct(ABxAC, xyz)<0 )
        {
            Membrane_triangle_list[i][1]=Point_C;
            Membrane_triangle_list[i][2]=Point_B;
            //cout<<"min"<<endl;
        }
    } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
}// END OF: Membrane_Normal_direction_Identifier function

//end of normal_direction_identifire***
