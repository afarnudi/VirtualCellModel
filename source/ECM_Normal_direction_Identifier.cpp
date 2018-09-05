//
//  ECM_Normal_direction_Identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.hpp"

void ECM::ECM_Normal_direction_Identifier(double x, double y, double z){
    
    double AC[3], AB[3], ABxAC[3], refernece[3];
    int Point_A, Point_B, Point_C;
    cout<<ECM_triangle_list.size()<<endl;
    for(  int i=0;i<ECM_triangle_list.size();i++)
    {
        Point_A=ECM_triangle_list[i][0];
        Point_B=ECM_triangle_list[i][1];
        Point_C=ECM_triangle_list[i][2];
        
        AB[0]=ECM_Node_Position[Point_B][0]-ECM_Node_Position[Point_A][0];
        AB[1]=ECM_Node_Position[Point_B][1]-ECM_Node_Position[Point_A][1];
        AB[2]=ECM_Node_Position[Point_B][2]-ECM_Node_Position[Point_A][2];
        
        AC[0]=ECM_Node_Position[Point_C][0]-ECM_Node_Position[Point_A][0];
        AC[1]=ECM_Node_Position[Point_C][1]-ECM_Node_Position[Point_A][1];
        AC[2]=ECM_Node_Position[Point_C][2]-ECM_Node_Position[Point_A][2];
        
        refernece[0]=ECM_Node_Position[Point_A][0]-x;
        refernece[1]=ECM_Node_Position[Point_A][1]-y;
        refernece[2]=ECM_Node_Position[Point_A][2]-z;
        
        crossvector(ABxAC, AB, AC);
        //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
        
        if(innerproduct(ABxAC, refernece)<0 )
        {
            ECM_triangle_list[i][1]=Point_C;
            ECM_triangle_list[i][2]=Point_B;
            //cout<<"min"<<endl;
        }
    } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
    }// END OF: Membrane_Normal_direction_Identifier function

//end of normal_direction_identifire***
