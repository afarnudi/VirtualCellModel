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
//                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "dihedral(p1,p2,p3,p4)"));
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
        } else if (dihedrals[i].type == potentialModelIndex.Model["DihedralArea"]) {
            auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
            if (DFs_item == DFs_classes.end()) {
                
                
                DFs_classes.insert(dihedrals[i].class_label);
                DFs_index++;
                
                string potential = "3*K_bend*distance(p2,p3)*(PI-dihedral(p1,p2,p3,p4))^2/(distance(p1,p2)*sin(angle(p1,p2,p3))+distance(p4,p2)*sin(angle(p4,p2,p3)))";
                
                DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, potential));
                DihedralForces[DFs_index]->addGlobalParameter("PI", M_PI);
                DihedralForces[DFs_index]->addPerBondParameter("K_bend");
//                DihedralForces[DFs_index]->addPerBondParameter("SponAngle");
                
                system.addForce(DihedralForces[DFs_index]);
                
            }
            
            
            vector<double> parameters;
            parameters.push_back(dihedrals[i].bendingStiffnessinKJ);
//            parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
            
            DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
        }
        
    }
}


string generate_Julicher1996_mean_curvature_potential(int order);

void set_mean_curvature_forces(MeanCurvature**                           mean_curvature_interactinos,
                               vector<OpenMM::CustomCompoundBondForce*> &MeanCurvatureForces,
                               OpenMM::System                           &system){
    
    int MCs_index = -1;
    for (int node_order=0; mean_curvature_interactinos[node_order][0].curvature_type != 2*EndOfList; ++node_order) {
//        cout<<"\n"<<node_order<<"\n";
        set <std::string> MCs_classes;
        
        for (int node_index=0; mean_curvature_interactinos[node_order][node_index].curvature_type != EndOfList; ++node_index) {
            
            if (mean_curvature_interactinos[node_order][node_index].curvature_type == potentialModelIndex.Model["Julicher1996"]) {
//                cout<<"h";
                string class_label =mean_curvature_interactinos[node_order][node_index].class_label+"O"+to_string(node_order);
                auto MCs_item = MCs_classes.find(class_label);
                if (MCs_item == MCs_classes.end()) {
                    
                    MCs_classes.insert(class_label);
                    MCs_index++;
                    
//                    cout<<node_order<<" "<<mean_curvature_interactinos[node_order][node_index].atoms.size()<<endl;
                    
                    string potential = generate_Julicher1996_mean_curvature_potential(node_order);
                    
                    MeanCurvatureForces.push_back(new OpenMM::CustomCompoundBondForce(node_order+1, potential));
                    
                    MeanCurvatureForces[MCs_index]->addGlobalParameter("PI", M_PI);
                    MeanCurvatureForces[MCs_index]->addPerBondParameter("k");
                    //                    MeanCurvatureForces[MCs_index]->addPerBondParameter("SponAngle");
                    
                    system.addForce(MeanCurvatureForces[MCs_index]);
                    
                }
                
                vector<double> parameters;
                parameters.push_back(mean_curvature_interactinos[node_order][node_index].curvatureStiffnessinKJpermol);
                //                parameters.push_back(dihedrals[i].spontaneousBendingAngleInRad);
//                cout<<MCs_index;
                MeanCurvatureForces[MCs_index]->addBond(mean_curvature_interactinos[node_order][node_index].atoms, parameters);
            }
        }
    }
    
}

string generate_Julicher1996_mean_curvature_potential(int node_order){
    string potential="0.75*k*(";
    string numerator="";
        cout<<node_order<<endl<<endl;
    
    
    
    numerator+="distance(p1,p"+to_string(2);
//    numerator+=")*dihedral(p"+to_string(node_order+1)+",p"+to_string(2)+",p1,p"+to_string(3)+")";
    numerator+=")*(PI-abs(dihedral(p"+to_string(node_order+1)+",p"+to_string(2)+",p1,p"+to_string(3)+")))";
    for (int i=3; i<node_order+1; i++) {
        numerator+="+distance(p1,p"+to_string(i);
//        numerator+=")*dihedral(p"+to_string(i-1)+",p"+to_string(i)+",p1,p"+to_string(i+1)+")";
        numerator+=")*(PI-abs(dihedral(p"+to_string(i-1)+",p"+to_string(i)+",p1,p"+to_string(i+1)+")))";
    }
    numerator+="+distance(p1,p"+to_string(node_order+1);
//    numerator+=")*dihedral(p"+to_string(node_order)+",p"+to_string(node_order+1)+",p1,p"+to_string(2)+")";
    numerator+=")*(PI-abs(dihedral(p"+to_string(node_order)+",p"+to_string(node_order+1)+",p1,p"+to_string(2)+")))";
    
    potential+=numerator+")^2/(";
    string denominator="";
    for (int i=2; i<node_order+1; i++) {
        denominator+="distance(p1,p"+to_string(i);
        denominator+=")*distance(p1,p"+to_string(i+1);
        denominator+=")*abs(sin(angle(p"+to_string(i)+",p1,p"+to_string(i+1)+")))+";
    }
    denominator+="distance(p1,p"+to_string(node_order+1);
    denominator+=")*distance(p1,p"+to_string(2);
    denominator+=")*abs(sin(angle(p"+to_string(node_order+1)+",p1,p"+to_string(2)+")))";
    
    potential+=denominator+")";
        cout<<potential<<endl<<endl;
    
    
    
    return potential;
}
