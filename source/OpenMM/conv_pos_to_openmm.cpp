#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_funcs.hpp"
#include "OpenMM_structs.h"
#include <string.h>

MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem) {
    const int mem_num_atom = mem.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[mem_num_atom];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int C=0;
    
    for (int i=0; i<mem_num_atom; i++) {
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
        myatominfo[i].radius=mem.get_node_radius();
        myatominfo[i].sigma_LJ_12_6=mem.get_sigma_LJ_12_6();
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
        myatominfo[i].radius=act.get_node_radius();
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
    
    for (int i=0; i<ecm_num_atom; i++) {
        myatominfo[i].type=C;
        myatominfo[i].class_label="ECM";
        std::string str = ecm.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].symbol = 'E';
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
        myatominfo[i].radius=ecm.get_node_radius();
        myatominfo[i].sigma_LJ_12_6=ecm.get_sigma_LJ_12_6();
        myatominfo[i].epsilon_LJ_12_6=ecm.get_epsilon_LJ_12_6();
        myatominfo[i].ext_force_model=ecm.get_ext_force_model();
        myatominfo[i].ext_force_constants[0]=ecm.get_kx();
        myatominfo[i].ext_force_constants[1]=ecm.get_ky();
        myatominfo[i].ext_force_constants[2]=ecm.get_kz();
        
    }
    //End of list
    //    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}

MyAtomInfo* convert_Chromatin_position_to_openmm(Chromatin chromo){
    const int chromo_num_atom = chromo.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[chromo_num_atom];
    
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
        
    }
    return myatominfo;
}
