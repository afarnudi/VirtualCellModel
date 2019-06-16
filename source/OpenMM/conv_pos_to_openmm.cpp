#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_funcs.hpp"
#include "OpenMM_structs.h"

MyAtomInfo* convert_membrane_position_to_openmm(Membrane mem) {
    const int mem_num_atom = mem.get_num_of_nodes();
    MyAtomInfo* myatominfo = new MyAtomInfo[mem_num_atom];
    
    //used in openmm to specify different types of atoms. I don't know what the application is at the moment.
    int C=0;
    
    for (int i=0; i<mem_num_atom; i++) {
        myatominfo[i].type=C;
        std::string str = mem.return_label();
        myatominfo[i].pdb = new char[str.length() + 1];
        strcpy(myatominfo[i].pdb, str.c_str());
        myatominfo[i].initPosInAng[0]=mem.get_node_position(i, 0);
        myatominfo[i].initPosInAng[1]=mem.get_node_position(i, 1);
        myatominfo[i].initPosInAng[2]=mem.get_node_position(i, 2);
        myatominfo[i].posInAng[0]=mem.get_node_position(i, 0);
        myatominfo[i].posInAng[1]=mem.get_node_position(i, 1);
        myatominfo[i].posInAng[2]=mem.get_node_position(i, 2);
        myatominfo[i].mass=mem.get_node_mass();
    }
    //End of list
//    myatominfo[mem_num_atom].type=-1;
    
    return myatominfo;
}
