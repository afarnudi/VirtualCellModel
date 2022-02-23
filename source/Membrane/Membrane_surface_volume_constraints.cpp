//
//  check.cpp
//  Mem
//
//  Created by Ali Farnudi on 24/09/2020.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Membrane.h"


void Membrane::assign_surface_volume_constraints(void){
    calculate_volume_and_surface_area();
    if (surface_constraint_model!=potentialModelIndex.Model["None"]) {
        if (SurfaceConstraintValue_stat=="Au") {
            if (surface_constraint_model==potentialModelIndex.Model["LocalConstraint"]) {
                SurfaceConstraintValue= surface_area/Num_of_Triangles;
            } else if(surface_constraint_model==potentialModelIndex.Model["GlobalConstraint"]) {
                SurfaceConstraintValue= surface_area;
            }
        } else{
            SurfaceConstraintValue= stod(SurfaceConstraintValue_stat);
        }
        if (LinearReducedSrfaceVolume==0) {
            SurfaceConstraintValue*=SurfaceConstraintRatio;
        }
        
    }
    if (volume_constraint_model!=potentialModelIndex.Model["None"]) {
        if (VolumeConstraintValue_stat != "Au") {
            VolumeConstraintValue = stod(VolumeConstraintValue_stat);
        } else {
            VolumeConstraintValue = volume*VolumeConstraintRatio;
            if (LinearReducedSrfaceVolume==0) {
                VolumeConstraintValue *= VolumeConstraintRatio;
            }
        }
    }
    
}
