#include "OpenMM_funcs.hpp"




Angles* convert_Chromatin_angle_bond_info_to_openmm(Chromatin chromo) {
    const int chromo_num_angles = chromo.get_num_of_angle_bonds();
    Angles* angles = new Angles[chromo_num_angles];
    
    
    bool COSpotential=false;
    
       
    for (int i=0; i<chromo_num_angles; i++) {
        angles[i].type = chromo.get_angle_spring_model();
        
        angles[i].atoms[0]=i;
        angles[i].atoms[1]=i+1;
        angles[i].atoms[2]=i+2;
        
        angles[i].class_label = chromo.get_label() + chromo.get_label();
        
        if (angles[i].type == potentialModelIndex.Model["angleCOS"]) {
            COSpotential=true;
            angles[i].spontaneousBendingAngleInRad = chromo.get_spontaneousBendingAngleInRad();
            angles[i].bendingStiffnessinKJpermol   = chromo.get_bendingStiffnessinKJpermol();
            
        }
        
    }
    if (COSpotential) {
        cout<<" cosine potential"<<endl;
        cout<<"\tBending rigidity (KJ.mol^-1) ="<<angles[0].bendingStiffnessinKJpermol<<endl;
        cout<<"\tSpontaneous bending angle (Rad) = "<<angles[0].bendingStiffnessinKJpermol<<endl;
        cout<<endl;
    } else {
        cout<<" None"<<endl;
    }
        
    return angles;
}
