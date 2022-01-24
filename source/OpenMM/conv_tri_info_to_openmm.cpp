#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem) {
    
    
    bool dihedralPotential = false;
    bool NonLinearDihedralPotential = false;
    bool meanCurvature = false;
    bool noPotential = false;
    
    const int mem_num_tri_pairs = mem.get_num_of_dihedral_elements();
    //    cout<<"**++**++**++\n\t"<<mem_num_tri_pairs<<"\n**++**++**++\n";
    Dihedrals* diatoms = new Dihedrals[mem_num_tri_pairs];
    
    double total_mem_area =  mem.get_surface_area();
    for (int i=0; i<mem_num_tri_pairs; i++) {
        //        vector<int> tri_pair_nodes(mem.get_traingle_pair_nodes_list(i));
        
        diatoms[i].type=mem.get_bending_model();
        if (diatoms[i].type == potentialModelIndex.Model["None"]){
            noPotential = true;
        } else  {
            vector<int> tri_pair_nodes(mem.get_dihedral_atoms_list(i));
            diatoms[i].atoms.resize(4);
            diatoms[i].atoms[0]=tri_pair_nodes[0];
            diatoms[i].atoms[1]=tri_pair_nodes[1];
            diatoms[i].atoms[2]=tri_pair_nodes[2];
            diatoms[i].atoms[3]=tri_pair_nodes[3];
            
            diatoms[i].bendingStiffnessinKJ = mem.get_bending_stiffness_coefficient();
            diatoms[i].class_label = mem.get_label();
            diatoms[i].spontaneousBendingAngleInRad = mem.get_spontaneous_angle_in_Rad(i);
            if (diatoms[i].type == potentialModelIndex.Model["Dihedral"]){
                dihedralPotential = true;
            } else if (diatoms[i].type == potentialModelIndex.Model["ExpDihedral"]){
                NonLinearDihedralPotential = true;
            }
//            else if (diatoms[i].type == potentialModelIndex.Model["cot_weight"]){
//                diatoms[i].total_mem_area = total_mem_area;
//                meanCurvature = true;
//            } else if (diatoms[i].type == potentialModelIndex.Model["DihedralArea"]){
////               diatoms[i].total_mem_area = total_mem_area;
//               meanCurvature = true;
//           }
        }
        
        
    }
    if (dihedralPotential) {
        cout<<" Cosine(dihedral)"<<endl;
        cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
    }
    if (NonLinearDihedralPotential) {
        cout<<" Non-linear exponential dihedral"<<endl;
        cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
    }
//    if (meanCurvature) {
//        cout<<" Mean curvature"<<endl;
//        cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
//    }
    
    if (noPotential || mem_num_tri_pairs==0) {
        cout<<" None"<<endl;
    }
    
    
    return diatoms;
}

Triangles* convert_membrane_triangle_info_to_openmm(Membrane &mem) {
    
    
    bool GlobalSurfacePotential = false;
    bool GlobalVolumePotential = false;
    bool LocalSurfacePotential = false;
    bool noSurfacePotential = false;
    bool noVolumePotential = false;
    
    const int mem_num_tris = mem.get_num_of_triangle();
    Triangles* triatoms = new Triangles[mem_num_tris];
    
    
    for (int i=0; i<mem_num_tris; i++) {
        
        triatoms[i].surface_type = mem.get_surface_constraint_model();
        triatoms[i].volume_type  = mem.get_volume_constraint_model();
        
        if (triatoms[i].surface_type != potentialModelIndex.Model["None"] || triatoms[i].volume_type != potentialModelIndex.Model["None"]) {
            
            vector<int> tri_nodes(mem.get_triangle_atoms_list(i));
            triatoms[i].atoms.resize(3);
            triatoms[i].atoms[0]=tri_nodes[0];
            triatoms[i].atoms[1]=tri_nodes[1];
            triatoms[i].atoms[2]=tri_nodes[2];
            
            if (triatoms[i].surface_type != potentialModelIndex.Model["None"]){
                if (triatoms[i].surface_type == potentialModelIndex.Model["LocalConstraint"]) {
                    LocalSurfacePotential = true;
                } else {
                    GlobalSurfacePotential = true;
                }
                triatoms[i].SurfaceConstraintStiffnessinKJpermolperNm2 = mem.get_surface_constraint_stiffness_coefficient();
                triatoms[i].class_label = mem.get_label();
                triatoms[i].SurfaceConstraintValue = mem.get_surface_constraint_area();
            } else {
                noSurfacePotential = true;
            }
            
            if (triatoms[i].volume_type != potentialModelIndex.Model["None"]){
                if (triatoms[i].volume_type == potentialModelIndex.Model["GlobalConstraint"]) {
                    GlobalVolumePotential = true;
                }
                triatoms[i].VolumeConstraintStiffnessinKJpermolperNm3 = mem.get_volume_constraint_stiffness_coefficient();
                triatoms[i].class_label = mem.get_label();
                triatoms[i].VolumeConstraintValue = mem.get_volume_constraint_value();
            } else {
                noVolumePotential = true;
            }
            
        } else {
            noSurfacePotential = true;
            noVolumePotential = true;
        }
        
        
    }
    cout<<"\tArea potential:";
    if (LocalSurfacePotential) {
        cout<<" Local surface constraint"<<endl;
        cout<<"\t\tCoeficient (KJ . mol^-1 . Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\t\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (GlobalSurfacePotential) {
        cout<<" Global surface constraint"<<endl;
        cout<<"\t\tCoeficient (KJ . mol^-1 . Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\t\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (noSurfacePotential || mem_num_tris==0) {
        cout<<" None"<<endl;
    }
    cout<<"\tVolume potential:";
    if (GlobalVolumePotential) {
        cout<<" Global volume constraint"<<endl;
        cout<<"\t\tCoeficient (KJ . mol^-1 .Nm^-3) = "<<mem.get_volume_constraint_stiffness_coefficient() <<endl;
        cout<<"\t\tVolume (Nm^3) = "<<mem.get_volume_constraint_value() <<endl;
    }
    if (noVolumePotential || mem_num_tris==0) {
        cout<<" None"<<endl;
    }
    
    //    exit(0);
    return triatoms;
}


MeanCurvature** convert_membrane_curvature_info_to_openmm(Membrane &mem) {
    
    
    bool curvaturePotential = false;
    bool noPotential = true;
    string node_orders;
    
//    const int mem_num_mcs;
    
    vector<vector<vector<int> > > nodeOrder_NodeIndex_NodeNeighbourList = mem.get_nodeOrder_NodeIndex_NodeNeighbourList();
    const int num_node_orders = nodeOrder_NodeIndex_NodeNeighbourList.size();
    
    MeanCurvature** mcatoms = new MeanCurvature*[num_node_orders];
    for (int node_order=0; node_order<num_node_orders; node_order++) {
        const int num_of_interactions =nodeOrder_NodeIndex_NodeNeighbourList[node_order].size();
        mcatoms[node_order] = new MeanCurvature[num_of_interactions];
        if (num_of_interactions!=0) {
            node_orders+=", "+to_string(node_order);
        }
        for (int node_index=0; node_index<num_of_interactions; node_index++) {
            mcatoms[node_order][node_index].curvature_type = mem.get_mean_curvature_model();
            if (mcatoms[node_order][node_index].curvature_type != potentialModelIndex.Model["None"] ) {
                if (mcatoms[node_order][node_index].curvature_type == potentialModelIndex.Model["Julicher1996"]) {
                    curvaturePotential = true;
                    noPotential = false;
                    
                    mcatoms[node_order][node_index].atoms=nodeOrder_NodeIndex_NodeNeighbourList[node_order][node_index];
                    mcatoms[node_order][node_index].curvatureStiffnessinKJpermol = mem.get_bending_stiffness_coefficient();
                    mcatoms[node_order][node_index].class_label = mem.get_label();
                }
                
                
            }
        }
        
        
        
        
        
    }
    
    if (curvaturePotential) {
        cout<<" Julicher (1996) descretisation"<<endl;
        cout<<"\tNode orders ="<<node_orders.erase(0,1)<<endl;
        cout<<"\tCoeficient (KJ / mol) = "<<mem.get_bending_stiffness_coefficient() <<endl;
//        cout<<"\tSpontaneous curvature = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (noPotential){
        cout<<" None"<<endl;
    }
    
    //    exit(0);
    return mcatoms;
}
