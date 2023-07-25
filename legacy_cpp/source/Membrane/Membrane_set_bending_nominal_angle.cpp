//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 24/09/2020.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::set_bending_nominal_angle(void){
    Triangle_pair_angles_in_radians.clear();
    Triangle_pair_angles_in_radians.resize(Num_of_Triangle_Pairs, M_PI);
    
    double avgangle=0;
    if (Triangle_pair_angle_stat=="Au" || Triangle_pair_angle_stat=="Av" ) {
        int pos1,pos2,pos3, pos4;
        double  N1[3], N2[3];
        double x1_3[3], x1_4[3], x2_3[3], x2_4[3];
        double  b1[3], b2[3], b3[3];
        //source: https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
        for (int i=0; i<Num_of_Triangle_Pairs; i++) {
            
            pos1=dihedral_atoms[i][0];
            pos2=dihedral_atoms[i][1];
            pos3=dihedral_atoms[i][2];
            pos4=dihedral_atoms[i][3];
            
            for (int index=0; index<3; index++) {
                b1[index]  = Node_Position[pos2][index]-Node_Position[pos1][index];
                b2[index]  = Node_Position[pos3][index]-Node_Position[pos2][index];
                b3[index]  = Node_Position[pos4][index]-Node_Position[pos3][index];
            }
            double b2length=vector_length(b2);
            double b1length=vector_length(b1);
            double b3length=vector_length(b3);
            
            for (int index=0; index<3; index++) {
                b2[index] /= b2length;
                b1[index] /= b1length;
                b3[index] /= b3length;
            }
            crossvector(N1, b1, b2);
            crossvector(N2, b2, b3);

            double x = innerproduct(N1, N2);
            double tempvec[3];
            crossvector(tempvec, N1, b2);
            double y = innerproduct(tempvec, N2);
            double angleInRad = atan2(y,x);
            angleInRad*=-1;
//            cout<<angleInRad*180/M_PI<<endl;
            
            if (Triangle_pair_angle_stat=="Au") {
                Triangle_pair_angles_in_radians[i]=angleInRad;
                
            } else {
                avgangle+=angleInRad;
            }
            
        }
        avgangle/=Num_of_Triangle_Pairs;
        if (Triangle_pair_angle_stat=="Au") {
            cout<<"Using mesh initial triangle pair angles as the nominal angle."<<endl;
        } else {
            for (int i=0; i<Num_of_Triangle_Pairs; i++) {
                Triangle_pair_angles_in_radians[i]=avgangle;
            }
            cout<<"Using the average mesh triangle pair angle ("<<avgangle<<" in radians) as the nominal angle."<<endl;
        }
        
    } else {
        cout<<"Using "<<Triangle_pair_Nominal_angle_in_degrees*M_PI/180.<<" (in radians) as the triangle pair nominal angle."<<endl;
        for (int i=0; i<Num_of_Triangle_Pairs; i++) {
            Triangle_pair_angles_in_radians[i] = Triangle_pair_Nominal_angle_in_degrees*M_PI/180.;
        }
    }

    
}
