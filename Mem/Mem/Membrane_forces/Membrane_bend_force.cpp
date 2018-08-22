//
//  Membrane_bend_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Membrane_bending_potetial(void){
    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    //    double temp_potential_energy = 0.0;
    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....
    // Beginning of the  triangle-triangle (bending) force calculations
    for(int i=0 ;i<Membrane_num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        pos1=Membrane_Triangle_Pair_Nodes[i][0];
        pos2=Membrane_Triangle_Pair_Nodes[i][3];
        pos3=Membrane_Triangle_Pair_Nodes[i][1];
        pos4=Membrane_Triangle_Pair_Nodes[i][2];
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            Ep2p1[index]=temp_p2[index]-temp_p1[index];
        }
        crossvector(N1, p3p1, p3p2);
        crossvector(N2, p4p2, p4p1);
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        crossvector(N3, N1, N2);
        double temp_Ep2p1_length=vectorlength(Ep2p1);
        //        if ((sinus- temp_sinus)>0.0001) {
        //            cout<< "oops!"<<endl;
        //            cout<<"temp_sinus="<<temp_sinus<<"\nsinus="<<sinus<<endl;
        //            exit (EXIT_FAILURE);
        //        }
        F0 = -(1-2*signbit(innerproduct(N3, Ep2p1)))*Membrane_bending_coefficient*sinus;
        //        cout<<"\nF0="<<F0<<endl;
        //        if( parallelORantiparallel(xpos ) == +1 )
        //        {
        //            F0=-F0;
        //        }
        double temp_N1_length_squared=vectorlength(N1)*vectorlength(N1);
        double temp_N2_length_squared=vectorlength(N2)*vectorlength(N2);
        //**************************************************** force calculation
//        F0*=-1;
        for (int l=0; l<3; l++) {
            F3[l]= F0 * temp_Ep2p1_length* N1[l]/ temp_N1_length_squared;
            
            F4[l]= F0 * temp_Ep2p1_length* N2[l]/ temp_N2_length_squared;
            
            F1[l]= (F0/temp_Ep2p1_length)*( innerproduct(p3p2,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p2,Ep2p1)*N2[l]/temp_N2_length_squared );
            
            F2[l]= (-F0/temp_Ep2p1_length)*( innerproduct(p3p1,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p1,Ep2p1)*N2[l]/temp_N2_length_squared );
            Membrane_Node_Force[pos1][l] += F1[l];
            Membrane_Node_Force[pos2][l] += F2[l];
            Membrane_Node_Force[pos3][l] += F3[l];
            Membrane_Node_Force[pos4][l] += F4[l];
            
        }
        
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: Not very sure about the calculations *****************************
        //***************** for 'Total_Potential_Energy. Need to double check *********************************
        //*******************************************************************************************************
        Total_Potential_Energy += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
    }  // end of 'for-i'
}
