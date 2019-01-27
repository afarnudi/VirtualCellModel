//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"

void Membrane::error_calculator_step_1(void){
    //If this method works we should move this loop into the MD_Evolution_2 to save on computational time.
    for (int i=0; i<Num_of_Nodes; i++) {
        Node_Position_lag_step[i][0]=Node_Position[i][0];
        Node_Position_lag_step[i][1]=Node_Position[i][1];
        Node_Position_lag_step[i][2]=Node_Position[i][2];
        
        Node_Force_lag_step[i][0]=Node_Force[i][0];
        Node_Force_lag_step[i][1]=Node_Force[i][1];
        Node_Force_lag_step[i][2]=Node_Force[i][2];
    }
}

void Membrane::error_calculator_step_2(void){
    double forc_2=0, momentum_2=0, force_derivative=0;
    double second_term=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        double a[3]={Node_Force[i][0], Node_Force[i][1], Node_Force[i][2]};
        double b[3]={Node_Velocity[i][0], Node_Velocity[i][1], Node_Velocity[i][2]};
        double c[3]={
            (Node_Force[i][0]-Node_Force_lag_step[i][0])/(Node_Position[i][0]-Node_Position_lag_step[i][0]),
            (Node_Force[i][1]-Node_Force_lag_step[i][1])/(Node_Position[i][1]-Node_Position_lag_step[i][1]),
            (Node_Force[i][2]-Node_Force_lag_step[i][2])/(Node_Position[i][2]-Node_Position_lag_step[i][2])};
        
        forc_2 += vector_length_squared(a);
        momentum_2 = vector_length_squared(b);
        force_derivative = vector_length(c);
        second_term+=momentum_2*force_derivative;
    }
    double dt=GenConst::MD_Time_Step;
    error=dt*dt*dt*(forc_2+second_term)/6.0;
//    error_com+=error;
}
