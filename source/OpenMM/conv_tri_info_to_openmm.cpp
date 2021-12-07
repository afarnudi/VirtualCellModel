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
                triatoms[i].VolumeConstraintValue = mem.get_volume_constraint_volume();
            } else {
                noVolumePotential = true;
            }
            
        } else {
            noSurfacePotential = true;
            noVolumePotential = true;
        }
        
        
    }
    cout<<"Area potential:";
    if (LocalSurfacePotential) {
        cout<<" Local surface constraint"<<endl;
        cout<<"\tCoeficient (KJ . mol . Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (GlobalSurfacePotential) {
        cout<<" Global surface constraint"<<endl;
        cout<<"\tCoeficient (KJ . mol . Nm^-2) = "<<mem.get_surface_constraint_stiffness_coefficient() <<endl;
        cout<<"\tSurface area (Nm^2) = "<<mem.get_surface_constraint_area() <<endl;
    }
    if (noSurfacePotential || mem_num_tris==0) {
        cout<<" None"<<endl;
    }
    cout<<"Volume potential:";
    if (GlobalVolumePotential) {
        cout<<" Global volume constraint"<<endl;
        cout<<"\tCoeficient (KJ . mol .Nm^-3) = "<<mem.get_volume_constraint_stiffness_coefficient() <<endl;
        cout<<"\tVolume (Nm^3) = "<<mem.get_volume_constraint_volume() <<endl;
    }
    if (noVolumePotential || mem_num_tris==0) {
        cout<<" None"<<endl;
    }
    
//    exit(0);
    return triatoms;
}
