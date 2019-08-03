#include "OpenMM_funcs.hpp"

Bonds* convert_Actin_bond_info_to_openmm(Actin act) {
    const int act_num_bonds = act.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[act_num_bonds];
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<act_num_bonds; i++) {
        bonds[i].type = act.get_spring_model();
        bonds[i].atoms[0]=act.get_node_pair(i, 0);
        bonds[i].atoms[1]=act.get_node_pair(i, 1);
        switch (bonds[i].type) {
                //FENE
            case 1:
//                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
//                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                bonds[i].nominalLengthInAngstroms=act.get_avg_node_dist();
                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                
                
        }
        
        
    }
    
    return bonds;
}

Bonds* convert_ECM_bond_info_to_openmm(ECM ecm) {
    const int ecm_num_bonds = ecm.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[ecm_num_bonds];
//    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
//    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<ecm_num_bonds; i++) {
        bonds[i].type = ecm.get_spring_model();
        bonds[i].atoms[0]=ecm.get_node_pair(i, 0);
        bonds[i].atoms[1]=ecm.get_node_pair(i, 1);
        switch (bonds[i].type) {
                //FENE
            case 1:
                ecm.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                bonds[i].stiffnessInKcalPerAngstrom2=ecm.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                bonds[i].nominalLengthInAngstroms=ecm.get_avg_node_dist();
                bonds[i].stiffnessInKcalPerAngstrom2=ecm.get_spring_stiffness_coefficient();
                break;
                
                
        }
        
        
    }

    return bonds;
}

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = mem.get_spring_model();
        bonds[i].atoms[0]=mem.get_node_pair(i, 0);
        bonds[i].atoms[1]=mem.get_node_pair(i, 1);
        bonds[i].class_label = mem.get_label();
        bonds[i].nominalLengthInAngstroms=mem.get_avg_node_dist();
        switch (bonds[i].type) {
            //FENE
            case 1:
                mem.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
            //Harmonic
            case 2:
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
            //HarmonicX4
            case 3:
                bonds[i].stiffnessInKcalPerAngstrom4=mem.get_spring_stiffness_coefficient();
                
            
        }
       
        
    }
    if(bonds[0].type==2){
    cout<<"bond potential: Harmonic "<<endl;
    cout<<"spring coeficient (KJ per Nanometer2) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }

    if (bonds[0].type == 1) {
        cout<< "bond potential: FENE"<<endl;
        cout<<"spring coeficient (KJ per Nanometer2) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
        cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;        
        cout<<"lmin(nm)\tlmax(nm)\tle0(nm)   \tle1(nm)\n"<<bonds[0].FENE_lmin * OpenMM::NmPerAngstrom<<"\t"<<bonds[0].FENE_lmax * OpenMM::NmPerAngstrom<<"\t"<<bonds[0].FENE_le0 * OpenMM::NmPerAngstrom<<"\t"<<bonds[0].FENE_le1 * OpenMM::NmPerAngstrom<<endl;
    }
    
    if(bonds[0].type == 3){
        cout<<"bond potential: X4Harmonic "<<endl;
        cout<<"spring coeficient (KJ per Nanometer4) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
        cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    cout<<endl;
    
    
    
    
    return bonds;
}
