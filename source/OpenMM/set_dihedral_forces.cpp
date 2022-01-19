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
        } else if (dihedrals[i].type == potentialModelIndex.Model["cot_weight"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
               
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
//                string potential = "K_bend*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle) )*0.5*abs(distance(p2,p3)*(distance(p1,p2)*sin(angle(p1,p2,p3))+distance(p2,p4)*sin(angle(p3,p2,p4))))/";
//                potential+=to_string(dihedrals[i].total_mem_area*3);

//                string potential = "0.5*abs(distance(p2,p3)*(distance(p1,p2)*sin(angle(p1,p2,p3))+distance(p2,p4)*sin(angle(p3,p2,p4))))/";
//                potential+=to_string(dihedrals[i].total_mem_area*3);
                
                string potential = "K_bend*(1-cos(dihedral(p1,p2,p3,p4)-SponAngle) )/(0.5*(cot(angle(p2,p1,p3))+cot(angle(p2,p4,p3))))";
                
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


string generate_itzykson1986_mean_curvature_potential();

void set_mean_curvature_forces(MeanCurvature**                           mean_curvature_interactinos,
                               vector<OpenMM::CustomCompoundBondForce*> &MeanCurvatureForces,
                               OpenMM::System                           &system){
    
    
    for (int node_order=0; mean_curvature_interactinos[node_order][0].curvature_type != 2*EndOfList; ++node_order) {
        set <std::string> MCs_classes;
        int MCs_index = -1;
        for (int node_index=0; mean_curvature_interactinos[node_order][node_index].curvature_type != EndOfList; ++node_index) {
            
            if (mean_curvature_interactinos[node_order][node_index].curvature_type == potentialModelIndex.Model["Julicher19960"]) {
                string class_label =mean_curvature_interactinos[node_order][node_index].class_label+"O"+to_string(node_order);
                auto MCs_item = MCs_classes.find(class_label);
                if (MCs_item == MCs_classes.end()) {
                    
                    MCs_classes.insert(class_label);
                    MCs_index++;
                    
                    string potential = generate_itzykson1986_mean_curvature_potential();
                    
                    MeanCurvatureForces.push_back(new OpenMM::CustomCompoundBondForce(node_order, potential));
                    
                    MeanCurvatureForces[MCs_index]->addPerBondParameter("K_b");
//                    MeanCurvatureForces[MCs_index]->addPerBondParameter("SponAngle");

                    system.addForce(MeanCurvatureForces[MCs_index]);
                    
                }
                
                vector<double> parameters;
                parameters.push_back(mean_curvature_interactinos[node_order][node_index].curvatureStiffnessinKJpermol);
//                parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
               
                MeanCurvatureForces[MCs_index]->addBond(mean_curvature_interactinos[node_order][node_index].atoms, parameters);
            }
        }
    }
    
}

string generate_itzykson1986_mean_curvature_potential(){
    string potential;
    return potential;
}
