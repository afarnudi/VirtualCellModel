//
//  Membrane_bend_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Bending_potetial_2(double sin_theta_0){
    int pos1,pos2,pos3,pos4;
    double  N1[3], N2[3], N3[3], sin, sign, F0;
    double E[3],x1_3[3], x1_4[3], x2_3[3], x2_4[3];
    
    for(int i=0 ;i<Num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        pos1=Triangle_Pair_Nodes[i][0];
        pos2=Triangle_Pair_Nodes[i][3];
        pos3=Triangle_Pair_Nodes[i][1];
        pos4=Triangle_Pair_Nodes[i][2];
        for (int index=0; index<3; index++) {
            E[index]=Node_Position[pos4][index]-Node_Position[pos3][index];
            x1_3[index]=Node_Position[pos1][index]-Node_Position[pos3][index];
            x1_4[index]=Node_Position[pos1][index]-Node_Position[pos4][index];
            x2_3[index]=Node_Position[pos2][index]-Node_Position[pos3][index];
            x2_4[index]=Node_Position[pos2][index]-Node_Position[pos4][index];
        }
        
        crossvector(N1, x1_3, x1_4);
        crossvector(N2, x2_4, x2_3);
        crossvector(N3, N1, N2);
        
        double E_length=vector_length(E), N1_length=vector_length(N1), N1_length_2=N1_length*N1_length, N2_length=vector_length(N2), N2_length_2=N2_length*N2_length;
        double N1dotN2=innerproduct(N1, N2)/(N1_length*N2_length);
        sin=(1.000001-N1dotN2)/2;
        if(sin<0){
            sin=-sin;
        }
        sin=sqrt(sin);
        sign=-(1-2*std::signbit(innerproduct(N3, E)));
        F0 = (sin-sin_theta_0)*sign*E_length*E_length*Bending_coefficient*sin/(N1_length+N2_length);
        
        double temp_1, temp_2, temp_3_1, temp_3_2, temp_4_1, temp_4_2;
        temp_1=E_length/N1_length_2;
        temp_2=E_length/N2_length_2;
        temp_3_1=innerproduct(x1_4, E)/(N1_length_2*E_length);
        temp_3_2=innerproduct(x2_4, E)/(N2_length_2*E_length);
        temp_4_1=innerproduct(x1_3, E)/(N1_length_2*E_length);
        temp_4_2=innerproduct(x2_3, E)/(N2_length_2*E_length);
        
        
        for (int l=0; l<3; l++) {
            Node_Force[pos1][l] -=  F0*temp_1*N1[l];
            Node_Force[pos2][l] -=  F0*temp_2*N2[l];
            Node_Force[pos3][l] -=  F0*temp_3_1*N1[l]+F0*temp_3_2*N2[l];
            Node_Force[pos4][l] -= -F0*temp_4_1*N1[l]-F0*temp_4_2*N2[l];
        }
    }  // end of 'for-i'
}

double Membrane::calculate_bending_energy(){
    int pos1,pos2,pos3,pos4;
    double  N1[3], N2[3];
    double x1_3[3], x1_4[3], x2_3[3], x2_4[3];
    
    Total_Bending_Energy=0;
    
    for(int i=0 ;i<Num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        pos1=Triangle_Pair_Nodes[i][0];
        pos2=Triangle_Pair_Nodes[i][3];
        pos3=Triangle_Pair_Nodes[i][1];
        pos4=Triangle_Pair_Nodes[i][2];
        for (int index=0; index<3; index++) {
            x1_3[index]=Node_Position[pos1][index]-Node_Position[pos3][index];
            x1_4[index]=Node_Position[pos1][index]-Node_Position[pos4][index];
            x2_3[index]=Node_Position[pos2][index]-Node_Position[pos3][index];
            x2_4[index]=Node_Position[pos2][index]-Node_Position[pos4][index];
        }
        
        crossvector(N1, x1_3, x1_4);
        crossvector(N2, x2_4, x2_3);
        
        double N1_length=vector_length(N1), N2_length=vector_length(N2);
        double N1dotN2 = innerproduct(N1, N2)/(N1_length*N2_length);
        Total_Bending_Energy += 1-N1dotN2;
    }  // end of 'for-i'
    Total_Bending_Energy *= Bending_coefficient;
    
    return Total_Bending_Energy;
}
