//
//  Membrane_spring_force.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Actin.h"


double Actin::Hookian (double distance, double initial_distance){
    Total_Potential_Energy += 0.5*Spring_coefficient*(distance-initial_distance)*(distance-initial_distance);
    return -Spring_coefficient*(distance-initial_distance);
}


double Actin::Kelvin(double distance, int bond_index){
    
    double initial_distance=Node_Bond_relaxed_length[bond_index];
    
    Total_Potential_Energy += 0.5*Spring_coefficient*(distance-initial_distance)*(distance-initial_distance);
    
    Node_Bond_relaxed_length[bond_index] +=  GenConst::Step_Size_In_Fs*( (distance-initial_distance)/distance )/Kelvin_Damping_Coefficient;
//    if (bond_index == 0) {
//        cout<<"Node_Bond_relaxed_length[bond_index]="<<Node_Bond_relaxed_length[bond_index]<<endl;
//    }
    
    return -Spring_coefficient*(distance-initial_distance);
}

double Actin::Maxwell(double distance, int bond_index){
    
    double initial_distance=Node_Bond_relaxed_length[bond_index];
    double gamma_0 = (distance/initial_distance) - 1;
    
    Node_Bond_relaxed_length[bond_index] +=  Dashpot_Viscosity*(distance-initial_distance)*(1-exp_tau);
    //    if (bond_index == 0) {
    //        cout<<"Node_Bond_relaxed_length[bond_index]="<<Node_Bond_relaxed_length[bond_index]<<endl;
    //    }
//    if (bond_index == 10) {
//        cout<< Node_Bond_relaxed_length[bond_index]<<endl;
//    }
    return -Spring_coefficient*gamma_0*exp_tau;
}


