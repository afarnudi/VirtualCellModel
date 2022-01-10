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
    //            DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(dihedral(p1,p2,p3,p4))"));
                
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
        } else if (dihedrals[i].type == potentialModelIndex.Model["meanCurvature"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                /**
                 Based on: EQ 15
                 Random Surface Discretizations and the Renormalization of the Bending Rigidity
                 G. Gompper et D.M. Kroll
                 J. Phys. I France, 6 10 (1996) 1305-1320
                 DOI: https://doi.org/10.1051/jp1:1996246
                */
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                string potential = "0.5*K_bend*( 4*( 1-cos(dihedral(p1,p2,p3,p4)-SponAngle) )/(cot(angle(p2,p1,p3))+cot(angle(p2,p4,p3))) )";
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, potential));
                
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
        } else if (dihedrals[i].type == potentialModelIndex.Model["meanCurvature2"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                /**
                 Based on: EQ 15
                 Random Surface Discretizations and the Renormalization of the Bending Rigidity
                 G. Gompper et D.M. Kroll
                 J. Phys. I France, 6 10 (1996) 1305-1320
                 DOI: https://doi.org/10.1051/jp1:1996246
                */
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                string potential = "2*K_bend*min( ( 1-cos(dihedral(p1,p2,p3,p4)-SponAngle) )/(cot(angle(p2,p1,p3))+cot(angle(p2,p4,p3))),( 1-cos(dihedral(p2,p4,p1,p3)-SponAngle) )/(cot(angle(p4,p2,p1))+cot(angle(p4,p3,p1))) )";
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, potential));
                
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
        }
        
    }
}
