#include "Membrane.h"

void Membrane::ConstantSurfaceForceLocalTriangles()
{
    
    int temp_node_A, temp_node_B, temp_node_C;
    double AB[3], AC[3], temp_AB_length_squared, temp_AC_length_squared, f0, temp_ABxAC[3];
    
    // cout <<surfacearea(x,tri )<<"   s0="<< s0<<endl;
    double s0_i=0.41*Node_radius*Node_radius; //1.732=3^0.5
    double s_i=0;
    
    
    for(  int i=0;i<Triangle_list.size();i++)
    {
        temp_node_A=Triangle_list[i][0];
        temp_node_B=Triangle_list[i][1];
        temp_node_C=Triangle_list[i][2];
        
        AB[0]=Node_Position[temp_node_B][0]-Node_Position[temp_node_A][0];
        AB[1]=Node_Position[temp_node_B][1]-Node_Position[temp_node_A][1];
        AB[2]=Node_Position[temp_node_B][2]-Node_Position[temp_node_A][2];
        
        AC[0]=Node_Position[temp_node_C][0]-Node_Position[temp_node_A][0];
        AC[1]=Node_Position[temp_node_C][1]-Node_Position[temp_node_A][1];
        AC[2]=Node_Position[temp_node_C][2]-Node_Position[temp_node_A][2];
        
        crossvector(temp_ABxAC, AB, AC);
        s_i = vectorlength(temp_ABxAC)/2.0;
        f0= K_surfaceConstant_local*(s_i -  s0_i )/2.0*vectorlength(temp_ABxAC);
        
        for (int j=0; j<3; j++) {
            temp_node_A=Triangle_list[i][j%3];
            temp_node_B=Triangle_list[i][(j+1)%3];
            temp_node_C=Triangle_list[i][(j+2)%3];
            
            AB[0]=Node_Position[temp_node_B][0]-Node_Position[temp_node_A][0];
            AB[1]=Node_Position[temp_node_B][1]-Node_Position[temp_node_A][1];
            AB[2]=Node_Position[temp_node_B][2]-Node_Position[temp_node_A][2];
            
            AC[0]=Node_Position[temp_node_C][0]-Node_Position[temp_node_A][0];
            AC[1]=Node_Position[temp_node_C][1]-Node_Position[temp_node_A][1];
            AC[2]=Node_Position[temp_node_C][2]-Node_Position[temp_node_A][2];
            
            temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
            temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
            
            Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
            Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
            Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
        }
    }
    
}