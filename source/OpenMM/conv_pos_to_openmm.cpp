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
        std::string str = mem.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].initPosInAng[0]=mem.get_node_position(i, 0);
        myatominfo[i].initPosInAng[1]=mem.get_node_position(i, 1);
        myatominfo[i].initPosInAng[2]=mem.get_node_position(i, 2);
        myatominfo[i].posInAng[0]=mem.get_node_position(i, 0);
        myatominfo[i].posInAng[1]=mem.get_node_position(i, 1);
        myatominfo[i].posInAng[2]=mem.get_node_position(i, 2);
        myatominfo[i].mass=mem.get_node_mass();
        myatominfo[i].radius=mem.get_node_radius();
        myatominfo[i].sigma_LJ_12_6=mem.get_sigma_LJ_12_6();
        myatominfo[i].epsilon_LJ_12_6=mem.get_epsilon_LJ_12_6();
        
    }
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
        std::string str = act.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].initPosInAng[0]=act.get_node_position(i, 0);
        myatominfo[i].initPosInAng[1]=act.get_node_position(i, 1);
        myatominfo[i].initPosInAng[2]=act.get_node_position(i, 2);
        //        cout<<act.get_node_position(i, 0)<<"\t"<<act.get_node_position(i, 1)<<"\t"<<act.get_node_position(i, 2)<<"\n";
        //        cout<<myatominfo[i].initPosInAng[0]<<"\t"<<myatominfo[i].initPosInAng[1]<<"\t"<<myatominfo[i].initPosInAng[2]<<"\n";
        myatominfo[i].posInAng[0]=act.get_node_position(i, 0);
        myatominfo[i].posInAng[1]=act.get_node_position(i, 1);
        myatominfo[i].posInAng[2]=act.get_node_position(i, 2);
        myatominfo[i].mass=act.get_node_mass();
        myatominfo[i].radius=act.get_node_radius();
        
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
        std::string str = ecm.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].initPosInAng[0]=ecm.get_node_position(i, 0);
        myatominfo[i].initPosInAng[1]=ecm.get_node_position(i, 1);
        myatominfo[i].initPosInAng[2]=ecm.get_node_position(i, 2);
        myatominfo[i].posInAng[0]=ecm.get_node_position(i, 0);
        myatominfo[i].posInAng[1]=ecm.get_node_position(i, 1);
        myatominfo[i].posInAng[2]=ecm.get_node_position(i, 2);
        myatominfo[i].mass=ecm.get_node_mass();
        myatominfo[i].radius=ecm.get_node_radius();
        myatominfo[i].sigma_LJ_12_6=ecm.get_sigma_LJ_12_6();
        myatominfo[i].epsilon_LJ_12_6=ecm.get_epsilon_LJ_12_6();
        
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
        std::string str = chromo.get_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].initPosInAng[0]=chromo.get_node_position(i, 0);
        myatominfo[i].initPosInAng[1]=chromo.get_node_position(i, 1);
        myatominfo[i].initPosInAng[2]=chromo.get_node_position(i, 2);
        myatominfo[i].posInAng[0]=chromo.get_node_position(i, 0);
        myatominfo[i].posInAng[1]=chromo.get_node_position(i, 1);
        myatominfo[i].posInAng[2]=chromo.get_node_position(i, 2);
        myatominfo[i].mass=chromo.get_node_mass();
        myatominfo[i].radius=chromo.get_node_radius();
        
    }
    //End of list
    //    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}
