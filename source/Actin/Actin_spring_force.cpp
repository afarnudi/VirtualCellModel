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
    Node_Bond_relaxed_length[bond_index] +=  GenConst::MD_Time_Step*( (distance-initial_distance)/distance )/Kelvin_Damping_Coefficient;
    
    return -Spring_coefficient*(distance-initial_distance);
}
