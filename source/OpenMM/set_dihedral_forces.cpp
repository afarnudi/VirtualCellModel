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
        
        auto DFs_item = DFs_classes.find(dihedrals[i].class_label);
        if (DFs_item == DFs_classes.end()) {
            
            DFs_classes.insert(dihedrals[i].class_label);
            DFs_index++;
            
            DihedralForces.push_back(new OpenMM::CustomCompoundBondForce(4, "K_bend*(1+cos(dihedral(p1,p2,p3,p4)))"));
            DihedralForces[DFs_index]->addPerBondParameter("K_bend");
            
            system.addForce(DihedralForces[DFs_index]);
            
        }
        
        vector<double> parameters={dihedrals[i].bendingStiffnessinKJ};
        DihedralForces[DFs_index]->addBond(dihedrals[i].atoms, parameters);
    }
}
