#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_funcs.hpp"
#include "OpenMM_structs.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>

MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem) {
    const int mem_num_atom = mem.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[mem_num_atom];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int C=0;
    int atom_count = 0;
    if (mem.LockOnSphere_stat || mem.LockOnEllipsoid_stat) {
        atom_count = 1;
        int last_atom_index = mem.get_num_of_nodes()-atom_count;
        myatominfo[last_atom_index].type=C;
        myatominfo[last_atom_index].class_label="Membrane";
        std::string str = mem.get_label();
        myatominfo[last_atom_index].pdb = new char[str.length() + 1];
        strcpy(myatominfo[last_atom_index].pdb, str.c_str());
        myatominfo[last_atom_index].energyInKJ = 0;
        myatominfo[last_atom_index].symbol = 'M';
        myatominfo[last_atom_index].initPosInNm[0]=mem.LockOnSphereCoordinate[0];
        myatominfo[last_atom_index].initPosInNm[1]=mem.LockOnSphereCoordinate[1];
        myatominfo[last_atom_index].initPosInNm[2]=mem.LockOnSphereCoordinate[2];
        myatominfo[last_atom_index].posInNm[0]=mem.LockOnSphereCoordinate[0];
        myatominfo[last_atom_index].posInNm[1]=mem.LockOnSphereCoordinate[1];
        myatominfo[last_atom_index].posInNm[2]=mem.LockOnSphereCoordinate[2];
        myatominfo[last_atom_index].velocityInNmperPs[0]=0;
        myatominfo[last_atom_index].velocityInNmperPs[1]=0;
        myatominfo[last_atom_index].velocityInNmperPs[2]=0;
        myatominfo[last_atom_index].mass = 0;
        myatominfo[last_atom_index].radius = 0;
        myatominfo[last_atom_index].sigma_LJ_12_6 = 0;
        myatominfo[last_atom_index].epsilon_LJ_12_6 = 0;
        myatominfo[last_atom_index].ext_force_model = 0;
        myatominfo[last_atom_index].ext_force_constants[0] = 0;
        myatominfo[last_atom_index].ext_force_constants[1] = 0;
        myatominfo[last_atom_index].ext_force_constants[2] = 0;
    }
    
    for (int i=0; i<mem_num_atom - atom_count; i++) {
        myatominfo[i].type=C;
        myatominfo[i].class_label="Membrane";
        std::string str = mem.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].energyInKJ = 0;
        myatominfo[i].symbol = 'M';
        myatominfo[i].initPosInNm[0]=mem.get_node_position(i, 0);
        myatominfo[i].initPosInNm[1]=mem.get_node_position(i, 1);
        myatominfo[i].initPosInNm[2]=mem.get_node_position(i, 2);
        myatominfo[i].posInNm[0]=mem.get_node_position(i, 0);
        myatominfo[i].posInNm[1]=mem.get_node_position(i, 1);
        myatominfo[i].posInNm[2]=mem.get_node_position(i, 2);
        myatominfo[i].velocityInNmperPs[0]=mem.get_node_velocity(i, 0);
        myatominfo[i].velocityInNmperPs[1]=mem.get_node_velocity(i, 1);
        myatominfo[i].velocityInNmperPs[2]=mem.get_node_velocity(i, 2);
        myatominfo[i].mass=mem.get_node_mass();
        myatominfo[i].radius=mem.get_node_radius(i);
        myatominfo[i].sigma_LJ_12_6=mem.get_node_radius(i);
        myatominfo[i].epsilon_LJ_12_6=mem.get_epsilon_LJ_12_6();
        myatominfo[i].ext_force_model=mem.get_ext_force_model();
        myatominfo[i].ext_force_constants[0]=mem.get_kx();
        myatominfo[i].ext_force_constants[1]=mem.get_ky();
        myatominfo[i].ext_force_constants[2]=mem.get_kz();
        
        
    }
   // myatominfo[0].mass=0;
   // myatominfo[1].mass=0;
    //myatominfo[2].mass=0;
    //myatominfo[3].mass=1;
    //myatominfo[2].mass=1;
    //End of list
//    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}

MyAtomInfo* convert_Actin_position_to_openmm(Actin act){
    const int act_num_atom = act.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[act_num_atom];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int C=0;
    
    for (int i=0; i<act_num_atom; i++) {
        myatominfo[i].type=C;
        myatominfo[i].class_label="Actin";
        std::string str = act.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].symbol = 'A';
        myatominfo[i].energyInKJ = 0;
        myatominfo[i].initPosInNm[0]=act.get_node_position(i, 0);
        myatominfo[i].initPosInNm[1]=act.get_node_position(i, 1);
        myatominfo[i].initPosInNm[2]=act.get_node_position(i, 2);
        myatominfo[i].velocityInNmperPs[0]=act.get_node_velocity(i, 0);
        myatominfo[i].velocityInNmperPs[1]=act.get_node_velocity(i, 1);
        myatominfo[i].velocityInNmperPs[2]=act.get_node_velocity(i, 2);
        myatominfo[i].posInNm[0]=act.get_node_position(i, 0);
        myatominfo[i].posInNm[1]=act.get_node_position(i, 1);
        myatominfo[i].posInNm[2]=act.get_node_position(i, 2);
        myatominfo[i].mass=act.get_node_mass();
        myatominfo[i].radius=act.get_node_radius(i);
        myatominfo[i].ext_force_model=act.get_ext_force_model();
        myatominfo[i].ext_force_constants[0]=act.get_kx();
        myatominfo[i].ext_force_constants[1]=act.get_ky();
        myatominfo[i].ext_force_constants[2]=act.get_kz();
        
    }
    //End of list
    //    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}

MyAtomInfo* convert_ECM_position_to_openmm(ECM ecm) {
    const int ecm_num_atom = ecm.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[ecm_num_atom];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int C=0;
    
    double min_x,min_y,min_z,max_x,max_z,max_y;
    min_x = 1000000;
    min_y = 1000000;
    min_z = 1000000;
    max_x = -1000000;
    max_z = -1000000;
    max_y = -1000000;
    
    srand(time(NULL));
    
    for (int i=0; i<ecm_num_atom; i++) {
        myatominfo[i].type=C;
        myatominfo[i].class_label="ECM";
        std::string str = ecm.get_label();
        //std::string str ("ec0a");
        std::string str2 ("ecm9");
        myatominfo[i].pdb = new char[str.length() + 1];
        //strcpy(myatominfo[i].pdb, str.c_str());
        //myatominfo[i].symbol = 'E';
        myatominfo[i].energyInKJ = 0;
        myatominfo[i].initPosInNm[0]=ecm.get_node_position(i, 0);
        myatominfo[i].initPosInNm[1]=ecm.get_node_position(i, 1);
        myatominfo[i].initPosInNm[2]=ecm.get_node_position(i, 2);
        myatominfo[i].posInNm[0]=ecm.get_node_position(i, 0);
        myatominfo[i].posInNm[1]=ecm.get_node_position(i, 1);
        myatominfo[i].posInNm[2]=ecm.get_node_position(i, 2);
        myatominfo[i].velocityInNmperPs[0]=ecm.get_node_velocity(i, 0);
        myatominfo[i].velocityInNmperPs[1]=ecm.get_node_velocity(i, 1);
        myatominfo[i].velocityInNmperPs[2]=ecm.get_node_velocity(i, 2);
        myatominfo[i].mass=ecm.get_node_mass();
        myatominfo[i].radius=ecm.get_node_radius(i);
        myatominfo[i].sigma_LJ_12_6=ecm.get_sigma_LJ_12_6();
        myatominfo[i].epsilon_LJ_12_6=ecm.get_epsilon_LJ_12_6();
        myatominfo[i].ext_force_model=ecm.get_ext_force_model();
        myatominfo[i].ext_force_constants[0]=ecm.get_kx();
        myatominfo[i].ext_force_constants[1]=ecm.get_ky();
        myatominfo[i].ext_force_constants[2]=ecm.get_kz();
        
        int type = ecm.get_receptor_type();
        double center_x = ecm.get_receptor_center_x();
        double center_y = ecm.get_receptor_center_y();
        double center_z = ecm.get_receptor_center_z();
        
        int r = rand() % 100 ;
        int density = int(100*(ecm.get_receptor_density() + (myatominfo[i].initPosInNm[0]-center_x) * ecm.get_receptor_gradient_x() + (myatominfo[i].initPosInNm[1]-center_y) * ecm.get_receptor_gradient_y()  + (myatominfo[i].initPosInNm[2]-center_z) * ecm.get_receptor_gradient_z() )) ;
        if(type==1)
        {
        if(r< density)
        {
            myatominfo[i].symbol = 'E';
            strcpy(myatominfo[i].pdb, str.c_str());
            
        }
        else{
            myatominfo[i].symbol = 'e';
            strcpy(myatominfo[i].pdb, str2.c_str());
        }
        }
        else if (type==2)
        {
            if( ((myatominfo[i].initPosInNm[0]-center_x) > 25) || ((myatominfo[i].initPosInNm[0]-center_x) < -25) )
            {
                myatominfo[i].symbol = 'E';
                strcpy(myatominfo[i].pdb, str.c_str());
            }
            else
            {
                myatominfo[i].symbol = 'e';
                strcpy(myatominfo[i].pdb, str2.c_str());
            }
        }
        
        
        if (myatominfo[i].initPosInNm[0]<min_x) {
            min_x = myatominfo[i].initPosInNm[0];
        }
        
        if (myatominfo[i].initPosInNm[1]<min_y) {
            min_y = myatominfo[i].initPosInNm[1];
        }
        
        if (myatominfo[i].initPosInNm[2]<min_z) {
            min_z = myatominfo[i].initPosInNm[2];
        }
        
        if (myatominfo[i].initPosInNm[0]>max_x) {
            max_x = myatominfo[i].initPosInNm[0];
        }
        
        if (myatominfo[i].initPosInNm[2]>max_z) {
            max_z = myatominfo[i].initPosInNm[2];
        }
        
        if (myatominfo[i].initPosInNm[1]>max_y) {
            max_y = myatominfo[i].initPosInNm[1];
        }
        
    }
    
    
    for (int i=0; i<ecm_num_atom; i++) {
       
        if(myatominfo[i].initPosInNm[0]< (min_x+0.00001))
        {
            myatominfo[i].mass = 0;
        }
        
        if(myatominfo[i].initPosInNm[2]< (min_z+0.00001))
        {
            myatominfo[i].mass = 0;
        }
        
        if(myatominfo[i].initPosInNm[1]< (min_y+0.00001))
        {
            if(min_y != max_y)
            {
            myatominfo[i].mass = 0;
            }
        }
        
        if(myatominfo[i].initPosInNm[0]> (max_x-0.00001))
        {
            myatominfo[i].mass = 0;
        }
        
        if(myatominfo[i].initPosInNm[2]> (max_z-0.00001))
        {
            myatominfo[i].mass = 0;
        }
        
        
        
    }
    
    
    
    
    //End of list
    //    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}

MyAtomInfo* convert_Chromatin_position_to_openmm(Chromatin chromo){
    const int chromo_num_atom = chromo.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[chromo_num_atom];
    
    int counter=0;
    vector<vector<int> > v_list;
    vector<vector<double> > v_weight_list;
    if (GenConst::ChromatinVirtualSites) {
        v_list = chromo.get_Vsite_and_bindings_list();
        v_weight_list = chromo.get_Vsite_binding_weight_list();
    }
    
    for (int i=0; i<chromo_num_atom; i++) {
        
        myatominfo[i].type=chromo.get_node_type(i);
        myatominfo[i].class_label="Chromatin";
        std::string str = chromo.get_label() + std::to_string(chromo.get_node_type(i));
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].energyInKJ = 0;
        myatominfo[i].symbol = 'C';
        myatominfo[i].sigma_LJ_12_6=chromo.get_sigma_LJ_12_6(chromo.get_node_type(i));
        myatominfo[i].epsilon_LJ_12_6=chromo.get_epsilon_LJ_12_6(chromo.get_node_type(i));
        myatominfo[i].initPosInNm[0]=chromo.get_node_position(i, 0);
        myatominfo[i].initPosInNm[1]=chromo.get_node_position(i, 1);
        myatominfo[i].initPosInNm[2]=chromo.get_node_position(i, 2);
        myatominfo[i].posInNm[0]=chromo.get_node_position(i, 0);
        myatominfo[i].posInNm[1]=chromo.get_node_position(i, 1);
        myatominfo[i].posInNm[2]=chromo.get_node_position(i, 2);
        myatominfo[i].velocityInNmperPs[0]=chromo.get_node_velocity(i, 0);
        myatominfo[i].velocityInNmperPs[1]=chromo.get_node_velocity(i, 1);
        myatominfo[i].velocityInNmperPs[2]=chromo.get_node_velocity(i, 2);
        myatominfo[i].mass=chromo.get_node_mass();
        myatominfo[i].radius=chromo.get_node_radius();
        
        if (GenConst::ChromatinVirtualSites) {
            if (counter < v_list.size()) {
                if (v_list[counter][0] == i) {
                    myatominfo[i].mass = 0;
                    
                    myatominfo[i].vsite_atoms[0] = v_list[counter][1];
                    myatominfo[i].vsite_atoms[1] = v_list[counter][2];
                    
                    myatominfo[i].Vsite_weights[0] = v_weight_list[counter][0];
                    myatominfo[i].Vsite_weights[1] = v_weight_list[counter][1];
                    
                    myatominfo[i].posInNm[0]=chromo.get_node_position(myatominfo[i].vsite_atoms[0], 0)*myatominfo[i].Vsite_weights[0] + chromo.get_node_position(myatominfo[i].vsite_atoms[1], 0)*myatominfo[i].Vsite_weights[1];
                    myatominfo[i].posInNm[1]=chromo.get_node_position(myatominfo[i].vsite_atoms[0], 1)*myatominfo[i].Vsite_weights[0] + chromo.get_node_position(myatominfo[i].vsite_atoms[1], 1)*myatominfo[i].Vsite_weights[1];
                    myatominfo[i].posInNm[2]=chromo.get_node_position(myatominfo[i].vsite_atoms[0], 2)*myatominfo[i].Vsite_weights[0] + chromo.get_node_position(myatominfo[i].vsite_atoms[1], 2)*myatominfo[i].Vsite_weights[1];
                    
                    myatominfo[i].initPosInNm[0]=myatominfo[i].posInNm[0];
                    myatominfo[i].initPosInNm[1]=myatominfo[i].posInNm[1];
                    myatominfo[i].initPosInNm[2]=myatominfo[i].posInNm[2];
                    
                    
                    
                    counter++;
                }
            }
            
            
        }
        
    }
    
    return myatominfo;
}
