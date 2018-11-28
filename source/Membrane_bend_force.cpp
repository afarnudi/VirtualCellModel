//
//  Membrane_bend_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::Bending_potetial(void){
    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    //    double temp_potential_energy = 0.0;
    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....
    // Beginning of the  triangle-triangle (bending) force calculations
    for(int i=0 ;i<Num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        pos1=Triangle_Pair_Nodes[i][0];
        pos2=Triangle_Pair_Nodes[i][3];
        pos3=Triangle_Pair_Nodes[i][1];
        pos4=Triangle_Pair_Nodes[i][2];
        for (int index=0; index<3; index++) {
            temp_p1[index]=Node_Position[pos1][index];
            temp_p2[index]=Node_Position[pos2][index];
            temp_p3[index]=Node_Position[pos3][index];
            temp_p4[index]=Node_Position[pos4][index];
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
        F0 = -(1-2*signbit(innerproduct(N3, Ep2p1)))*Bending_coefficient*sinus;
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
            Node_Force[pos1][l] += F1[l];
            Node_Force[pos2][l] += F2[l];
            Node_Force[pos3][l] += F3[l];
            Node_Force[pos4][l] += F4[l];
            
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
        Total_Potential_Energy += Bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
    }  // end of 'for-i'
}

void Membrane::Bending_potetial_2(double sin_theta_0){
    int pos1,pos2,pos3,pos4;
    double  N1[3], N2[3], N3[3], sinus, sign, F0;
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
        
        double E_length=vectorlength(E), N1_length=vectorlength(N1), N1_length_2=N1_length*N1_length, N2_length=vectorlength(N2), N2_length_2=N2_length*N2_length;
        double N1dotN2=innerproduct(N1, N2)/(N1_length*N2_length);
        sinus=(1.000001-N1dotN2)/2;
		if(sinus<0){
			sinus=-sinus;
			}
		sinus=sqrt(sinus);
        sign=-(1-2*signbit(innerproduct(N3, E)));
        F0 = -(sinus-sin_theta_0)*sign*E_length*E_length*Bending_coefficient*sinus/(N1_length+N2_length);
        
        double temp_1, temp_2, temp_3_1, temp_3_2, temp_4_1, temp_4_2;
        temp_1=E_length/N1_length_2;
        temp_2=E_length/N2_length_2;
        temp_3_1=innerproduct(x1_4, E)/(N1_length_2*E_length);
        temp_3_2=innerproduct(x2_4, E)/(N2_length_2*E_length);
        temp_4_1=innerproduct(x1_3, E)/(N1_length_2*E_length);
        temp_4_2=innerproduct(x2_3, E)/(N2_length_2*E_length);
        
        
        for (int l=0; l<3; l++) {
            Node_Force[pos1][l] +=  F0*temp_1*N1[l];
            Node_Force[pos2][l] +=  F0*temp_2*N2[l];
            Node_Force[pos3][l] +=  F0*temp_3_1*N1[l]+F0*temp_3_2*N2[l];
            Node_Force[pos4][l] += -F0*temp_4_1*N1[l]-F0*temp_4_2*N2[l];
            
        }
    }  // end of 'for-i'
}
