#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

void set_angle_forces_minimum(Angles*                             angles,
                              vector<OpenMM::CustomAngleForce*>  &AngleForces,
                              OpenMM::System                     &system
                              ){
    
    for (int i=0; angles[i].type != EndOfList; ++i) {
//        cout<<i<<" "<<flush;
        
        if (angles[i].type == potentialModelIndex.Model["angleCOS"]) {
            
            if (AngleForces.size() == 0) {
                string k_bend = "K_bend";
                string theta0 = "theta0Bend";
                AngleForces.push_back(new OpenMM::CustomAngleForce(k_bend+"*(1-cos(theta-"+theta0+") )"));
                
//                AngleForces[0]->addPerAngleParameter(k_bend);
//                AngleForces[0]->addPerAngleParameter(theta0);
                AngleForces[0]->addGlobalParameter(k_bend, angles[i].bendingStiffnessinKJpermol);
                AngleForces[0]->addGlobalParameter(theta0, angles[i].spontaneousBendingAngleInRad);
                
                system.addForce(AngleForces[0]);
                
            }
//            vector<double> parameters={angles[i].bendingStiffnessinKJpermol,
//                                       angles[i].spontaneousBendingAngleInRad
//                                       };
//            AngleForces[0]->addAngle(angles[i].atoms[0],angles[i].atoms[1],angles[i].atoms[2],parameters);
            
            AngleForces[0]->addAngle(angles[i].atoms[0],angles[i].atoms[1],angles[i].atoms[2]);
        }
        
    }
}



void set_angle_forces(Angles*                             angles,
                      vector<OpenMM::CustomAngleForce*>  &AngleForces,
                      OpenMM::System                     &system
                      ){
    //AFs = AngleForces
    set <std::string> AFs_classes;
    int AFs_index = -1;
    
    for (int i=0; angles[i].type != EndOfList; ++i) {
//        cout<<i<<" "<<flush;
        
        if (angles[i].type == potentialModelIndex.Model["angleCOS"]) {
            auto AFs_item = AFs_classes.find(angles[i].class_label);
            if (AFs_item == AFs_classes.end()) {
                
                AFs_classes.insert(angles[i].class_label);
                AFs_index++;
                
                AngleForces.push_back(new OpenMM::CustomAngleForce("K_bend"+angles[i].class_label+"*(1-cos(theta-theta0"+angles[i].class_label+") )"));
    //            DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(dihedral(p1,p2,p3,p4))"));
                
                AngleForces[AFs_index]->addGlobalParameter("K_bend"+angles[i].class_label, angles[i].bendingStiffnessinKJpermol);
                AngleForces[AFs_index]->addGlobalParameter("theta0"+angles[i].class_label, angles[i].spontaneousBendingAngleInRad);
                
                
                system.addForce(AngleForces[AFs_index]);
                
            }
           
            AngleForces[AFs_index]->addAngle(angles[i].atoms[0],angles[i].atoms[1],angles[i].atoms[2]);
        }
        
    }
}
