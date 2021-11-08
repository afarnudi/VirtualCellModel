#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"

Dihedrals* convert_membrane_dihedral_info_to_openmm(Membrane &mem) {
    
    
    bool dihedralPotential = false;
    bool noPotential = false;
    
    const int mem_num_tri_pairs = mem.get_num_of_dihedral_elements();
//    cout<<"**++**++**++\n\t"<<mem_num_tri_pairs<<"\n**++**++**++\n";
    Dihedrals* diatoms = new Dihedrals[mem_num_tri_pairs];
    
    
    for (int i=0; i<mem_num_tri_pairs; i++) {
//        vector<int> tri_pair_nodes(mem.get_traingle_pair_nodes_list(i));
        
        diatoms[i].type=mem.get_bending_model();
        if (diatoms[i].type == potentialModelIndex.Model["Dihedral"]) {
            dihedralPotential = true;
            vector<int> tri_pair_nodes(mem.get_dihedral_atoms_list(i));
            diatoms[i].atoms.resize(4);
            diatoms[i].atoms[0]=tri_pair_nodes[0];
            diatoms[i].atoms[1]=tri_pair_nodes[1];
            diatoms[i].atoms[2]=tri_pair_nodes[2];
            diatoms[i].atoms[3]=tri_pair_nodes[3];

            diatoms[i].bendingStiffnessinKJ = mem.get_bending_stiffness_coefficient();
            diatoms[i].class_label = mem.get_label();
            diatoms[i].spontaneousBendingAngleInRad = mem.get_spontaneous_angle_in_Rad(i);
        } else if (diatoms[i].type == potentialModelIndex.Model["None"]){
            noPotential = true;
        }
        
        
    }
    if (dihedralPotential) {
        cout<<" Cosine(dihedral)"<<endl;
        cout<<"\tCoeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
    }
    
    if (noPotential || mem_num_tri_pairs==0) {
        cout<<" None"<<endl;
    }
    
    
    return diatoms;
}

Triangles* convert_membrane_triangle_info_to_openmm(Membrane &mem) {
    
    
    bool GlobalSurfacePotential = false;
    bool LocalSurfacePotential = false;
    bool noPotential = false;
    
    const int mem_num_tris = mem.get_num_of_triangle();
//    cout<<"**++**++**++\n\t"<<mem_num_tri_pairs<<"\n**++**++**++\n";
    Triangles* triatoms = new Triangles[mem_num_tris];
    
    
    for (int i=0; i<mem_num_tris; i++) {
//        vector<int> tri_pair_nodes(mem.get_traingle_pair_nodes_list(i));
        
        triatoms[i].type = mem.get_surface_constraint_model();
        if (triatoms[i].type == potentialModelIndex.Model["LocalConstraint"] || triatoms[i].type == potentialModelIndex.Model["GlobalConstraint"]) {
            if (triatoms[i].type == potentialModelIndex.Model["LocalConstraint"]) {
                LocalSurfacePotential = true;
            } else {
                GlobalSurfacePotential = true;
            }
            vector<int> tri_nodes(mem.get_triangle_atoms_list(i));
            triatoms[i].atoms.resize(3);
            triatoms[i].atoms[0]=tri_nodes[0];
            triatoms[i].atoms[1]=tri_nodes[1];
            triatoms[i].atoms[2]=tri_nodes[2];
            
            triatoms[i].ConstraintStiffnessinKJpermolperNm2 = mem.get_surface_constraint_stiffness_coefficient();
            triatoms[i].class_label = mem.get_label();
            triatoms[i].SurfaceConstraintValue = mem.get_surface_constraint_area();
        } else if (triatoms[i].type == potentialModelIndex.Model["None"]){
            noPotential = true;
        }
        
        
    }
    
    if (LocalSurfacePotential) {
        cout<<" Local surface constraint"<<endl;
        cout<<"\tCoeficient (KJ . mol^-2 .Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (GlobalSurfacePotential) {
        cout<<" Global surface constraint"<<endl;
        cout<<"\tCoeficient (KJ . mol^-2 .Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    
    if (noPotential || mem_num_tris==0) {
        cout<<" None"<<endl;
    }
    
    
    return triatoms;
}
