#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

void set_dihedral_forces(Dihedrals*                                 dihedrals,
                         vector<OpenMM::CustomCompoundBondForce*>  &DihedralForces,
                         OpenMM::System                            &system
                         ){
    //DFs = DihedralForces
    set <std::string> DFs_classes;
    int DFs_index = -1;
    
    for (int i=0; dihedrals[i].type != EndOfList; ++i) {
        
        if (dihedrals[i].type == potentialModelIndex.Model["Dihedral"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle) )"));
                
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                //
                //                if (GenConst::Periodic_box) {
                //                    DihedralForces[DFs_index]->setUsesPeriodicBoundaryConditions(true);
                //                }
                
                system.addForce(DihedralForces[DFs_index]);
            }
            
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        } else if (dihedrals[i].type == potentialModelIndex.Model["ExpDihedral"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "0.5*K_bend*(exp(2*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle))) -1)"));
                
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                //
                //                if (GenConst::Periodic_box) {
                //                    DihedralForces[DFs_index]->setUsesPeriodicBoundaryConditions(true);
                //                }
                
                system.addForce(DihedralForces[DFs_index]);
                
            }
            
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        } else if (dihedrals[i].type == potentialModelIndex.Model["SmoothEXP"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "0.04*K_bend*(exp(10*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle))) -1)"));
                
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                //
                //                if (GenConst::Periodic_box) {
                //                    DihedralForces[DFs_index]->setUsesPeriodicBoundaryConditions(true);
                //                }
                
                system.addForce(DihedralForces[DFs_index]);
                if (generalParameters.WantForce && generalParameters.force_group_count<31) {
                    generalParameters.force_group_count++;
                    DihedralForces[DFs_index]->setForceGroup(generalParameters.force_group_count);
                    string label = dihedrals[i].class_label+"_ExpDihedral_"+ to_string(generalParameters.force_group_count);
                    generalParameters.force_group_label.push_back(label);
                    
                }
            }
            
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        } else if (dihedrals[i].type == potentialModelIndex.Model["SmoothTheta4"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle))^2"));
                
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                //
                //                if (GenConst::Periodic_box) {
                //                    DihedralForces[DFs_index]->setUsesPeriodicBoundaryConditions(true);
                //                }
                
                system.addForce(DihedralForces[DFs_index]);
                if (generalParameters.WantForce && generalParameters.force_group_count<31) {
                    generalParameters.force_group_count++;
                    DihedralForces[DFs_index]->setForceGroup(generalParameters.force_group_count);
                    string label = dihedrals[i].class_label+"_Theta4_"+ to_string(generalParameters.force_group_count);
                    generalParameters.force_group_label.push_back(label);
                    
                }
                
            }
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        } else if (dihedrals[i].type == potentialModelIndex.Model["SmoothEXP46"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            
            if (DFs_item == DFs_classes.end()) {
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "0.5*K_bend*((exp(2*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle))) -1)-(dihedral(p1,p2,p3,p4)-SponAngle)^2)"));
                
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                //
                //                if (GenConst::Periodic_box) {
                //                    DihedralForces[DFs_index]->setUsesPeriodicBoundaryConditions(true);
                //                }
                
                system.addForce(DihedralForces[DFs_index]);
                if (generalParameters.WantForce && generalParameters.force_group_count<31) {
                    generalParameters.force_group_count++;
                    DihedralForces[DFs_index]->setForceGroup(generalParameters.force_group_count);
                    string label = dihedrals[i].class_label+"_Exp46_"+ to_string(generalParameters.force_group_count);
                    generalParameters.force_group_label.push_back(label);
                    
                }
                
            }
            
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        }
        
        
        
    }
}

