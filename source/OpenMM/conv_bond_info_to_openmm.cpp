#include "OpenMM_funcs.hpp"

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = mem.get_spring_model();
        bonds[i].atoms[0]=mem.get_node_pair(i, 0);
        bonds[i].atoms[1]=mem.get_node_pair(i, 1);
        bonds[i].class_label = mem.get_label() + mem.get_label();
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
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
                //Kelvin-Voigt
            case 4:
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
                
            case 5:
                bonds[i].stiffnessInKcalPerAngstrom2=mem.get_spring_stiffness_coefficient();
                break;
        }
        
        
    }
    if(bonds[0].type==2){
        cout<<"bond potential: Harmonic "<<endl;
        cout<<"spring coeficient (KJ per Nanometer2) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
        cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal <<endl;
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
    
    if(bonds[0].type == 4){
        cout<<"Membrane bond potential: Kelvin-Voigt "<<endl;
        cout<<"spring coeficient (KJ per Nanometer2) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm <<endl;
        cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    if(bonds[0].type==5){
        cout<<"bond potential: realHarmonic "<<endl;
        cout<<"spring coeficient (KJ per Nanometer2) ="<<mem.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
        cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    
    cout<<endl;
    
    
    
    
    return bonds;
}

Bonds* convert_Actin_bond_info_to_openmm(Actin act,MyAtomInfo* atoms) {
    const int act_num_bonds = act.get_num_of_node_pairs();
    const int act_abp_bonds = act.get_num_of_abp_pairs();
    const int act_MT_bonds = act.get_num_of_MT_pairs();
    Bonds* bonds = new Bonds[4*act_num_bonds + 4*act_abp_bonds + 4*act_MT_bonds];
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<act_num_bonds; i++) {
        bonds[i].type = act.get_spring_model();
        bonds[i].atoms[0]=act.get_node_pair(i, 0);
        bonds[i].atoms[1]=act.get_node_pair(i, 1);
        bonds[i].class_label = act.get_label() + act.get_label();
        switch (bonds[i].type) {
                //FENE
            case 1:
                //                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                //                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                //bonds[i].nominalLengthInAngstroms=act.get_avg_node_dist();
                bonds[i].nominalLengthInAngstroms=act.get_act_relaxlength(i) * act.get_act_r0factor();
                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                
                //Kelvin-Voigt
            case 4:
                //bonds[i].nominalLengthInAngstroms=act.get_avg_node_dist();
                bonds[i].nominalLengthInAngstroms=act.get_act_relaxlength(i) * act.get_act_r0factor();
                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                bonds[i].dampInKcalPsPerAngstrom2=act.get_kelvin_damping_coefficient();
                
//                if( (0.5*(atoms[bonds[i].atoms[0]].initPosInAng[0] + atoms[bonds[i].atoms[1]].initPosInAng[0]) ) >0  )
//                {
//                    bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
//
//                }
//                else
//                {
//                    bonds[i].stiffnessInKcalPerAngstrom2=2*act.get_spring_stiffness_coefficient();
//
//                }
                
                
                break;
                
        }
        
        
    }
    
    
    
    
    //Contractile element parallel to bond spring
    for (int i=act_num_bonds; i<(2*act_num_bonds); i++) {
        int model = act.get_contractile_model();
        switch (model) {
            case 1:
                bonds[i].type = 6;
                break;
            case 2:
                bonds[i].type = 8;
                break;
            case 3:
                bonds[i].type = 9;
                break;
        }
        bonds[i].k_F0 = 4;
        //bonds[i].type = 6;
        bonds[i].atoms[0]=act.get_node_pair(i-act_num_bonds, 0);
        bonds[i].atoms[1]=act.get_node_pair(i-act_num_bonds, 1);
        bonds[i].class_label = act.get_label() + act.get_label();
        
        //bonds[i].nominalLengthInAngstroms=act.get_avg_node_dist();
        bonds[i].nominalLengthInAngstroms=act.get_act_relaxlength(i-act_num_bonds) * act.get_act_r0factor();
        bonds[i].F0 = act.get_contractile_force();
        bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInAngstroms ;
        bonds[i].r_max = act.get_contractile_rmax() * bonds[i].nominalLengthInAngstroms;
        //bonds[i].hill_co = act.get_hill_co();
        bonds[i].hill_co = act.get_contractile_hill_co();
    }
    
    //Contractile spring
    for (int i=2*act_num_bonds; i<(3*act_num_bonds); i++) {
        bonds[i].type = 7;
        bonds[i+act_num_bonds].type = 7;
        bonds[i].atoms[0]=act.get_node_pair(i-2*act_num_bonds, 0);
        bonds[i].atoms[1]=act.get_node_pair(i-2*act_num_bonds, 1);
        bonds[i+act_num_bonds].atoms[0]=act.get_node_pair(i-2*act_num_bonds, 0);
        bonds[i+act_num_bonds].atoms[1]=act.get_node_pair(i-2*act_num_bonds, 1);
        bonds[i].class_label = act.get_label() + act.get_label();
        bonds[i+act_num_bonds].class_label = act.get_label() + act.get_label();
        
        //bonds[i].nominalLengthInAngstroms=act.get_avg_node_dist();
        bonds[i].nominalLengthInAngstroms=act.get_act_relaxlength(i-2*act_num_bonds) * act.get_act_r0factor();
        bonds[i].stiffnessInKcalPerAngstrom2 = act.get_contractile_k1();
        bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInAngstroms;
        bonds[i].r_max = bonds[i].nominalLengthInAngstroms;
        //bonds[i].r_max = act.get_avg_node_dist();
        
        bonds[i+act_num_bonds].nominalLengthInAngstroms=act.get_act_relaxlength(i-2*act_num_bonds) * act.get_act_r0factor();
        bonds[i+act_num_bonds].stiffnessInKcalPerAngstrom2 = act.get_contractile_k2();
        bonds[i+act_num_bonds].r_min = bonds[i+act_num_bonds].nominalLengthInAngstroms;
        bonds[i+act_num_bonds].r_max = act.get_contractile_rmax() * bonds[i+act_num_bonds].nominalLengthInAngstroms;
    }
    
    
    int n1=0;
    int n2=0;
    
    for (int i=4*act_num_bonds; i<(4*act_num_bonds + act_abp_bonds); i++) {
       
        int model = act.get_abp_model();
        switch (model) {
            case 1:
                bonds[i].type = 6;
                break;
            case 2:
                bonds[i].type = 8;
                break;
            case 3:
                bonds[i].type = 9;
                break;
        }
        
        bonds[i].k_F0 = 4;

        
        //bonds[i].type = 6;
        bonds[i].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i].class_label = act.get_label() + act.get_label();
        
        bonds[i].nominalLengthInAngstroms=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        
        
        
        bonds[i].F0 = act.get_abp_force();
        
//        if( (0.5*(atoms[bonds[i].atoms[0]].initPosInAng[0] + atoms[bonds[i].atoms[1]].initPosInAng[0]) ) >0  )
//        {
//            bonds[i].F0 =  act.get_abp_force();
//            n1++;
//        }
//        else
//        {
//            bonds[i].F0 = 10*act.get_abp_force();
//            n2++;
//        }
        
        
        bonds[i].r_min = act.get_abp_rmin() * bonds[i].nominalLengthInAngstroms;
        bonds[i].r_max = act.get_abp_rmax() * bonds[i].nominalLengthInAngstroms;
        bonds[i].hill_co = act.get_abp_hill_co();
        
        
        
        bonds[i+act_abp_bonds].type = 7;
        bonds[i+act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i+act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i+act_abp_bonds].class_label = act.get_label() + act.get_label();
        
        bonds[i+act_abp_bonds].nominalLengthInAngstroms=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        bonds[i+act_abp_bonds].stiffnessInKcalPerAngstrom2 = act.get_abp_k1();
        bonds[i+act_abp_bonds].r_min = act.get_abp_rmin() * bonds[i+act_abp_bonds].nominalLengthInAngstroms;
        bonds[i+act_abp_bonds].r_max = bonds[i+act_abp_bonds].nominalLengthInAngstroms;
        
        
        
        bonds[i+2*act_abp_bonds].type = 7;
        bonds[i+2*act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i+2*act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i+2*act_abp_bonds].class_label = act.get_label() + act.get_label();
        
        bonds[i+2*act_abp_bonds].nominalLengthInAngstroms=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        bonds[i+2*act_abp_bonds].stiffnessInKcalPerAngstrom2 = act.get_abp_k2();
        bonds[i+2*act_abp_bonds].r_min = bonds[i+2*act_abp_bonds].nominalLengthInAngstroms;
        bonds[i+2*act_abp_bonds].r_max = act.get_abp_rmax() * bonds[i+2*act_abp_bonds].nominalLengthInAngstroms;
        
        
        
        
        
        bonds[i+3*act_abp_bonds].type = act.get_abp_spring_model();
        bonds[i+3*act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i+3*act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i+3*act_abp_bonds].class_label = act.get_label() + act.get_label();
        switch (bonds[i+3*act_abp_bonds].type) {
                //FENE
            case 1:
                //                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                //                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                bonds[i+3*act_abp_bonds].nominalLengthInAngstroms=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKcalPerAngstrom2=act.get_abp_spring_stiffness_coefficient();
                break;
                
                //Kelvin-Voigt
            case 4:
                bonds[i+3*act_abp_bonds].nominalLengthInAngstroms=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKcalPerAngstrom2=act.get_abp_spring_stiffness_coefficient();
                bonds[i+3*act_abp_bonds].dampInKcalPsPerAngstrom2=act.get_kelvin_damping_coefficient();
                break;
                
        }
        
        
    }
    
    
    
    
    
    
    
    for (int i=4*act_num_bonds+ 4*act_abp_bonds; i<(4*act_num_bonds + 4*act_abp_bonds + act_MT_bonds); i++) {
           
            int model = act.get_MT_model();
            switch (model) {
                case 1:
                    bonds[i].type = 6;
                    break;
                case 2:
                    bonds[i].type = 8;
                    break;
                case 3:
                    bonds[i].type = 9;
                    break;
            }
            
            bonds[i].k_F0 = 4;

            
            //bonds[i].type = 6;
            bonds[i].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i].class_label = act.get_label() + act.get_label();
            
            bonds[i].nominalLengthInAngstroms=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            
            
            
           // bonds[i].F0 = act.get_MT_force();
            
            if( (0.5*(atoms[bonds[i].atoms[0]].initPosInAng[0] + atoms[bonds[i].atoms[1]].initPosInAng[0]) ) >0  )
            {
                bonds[i].F0 =  act.get_MT_force();
                n1++;
            }
            else
            {
                bonds[i].F0 = 5*act.get_MT_force();
                n2++;
            }
            
            
            bonds[i].r_min = act.get_MT_rmin() * bonds[i].nominalLengthInAngstroms;
            bonds[i].r_max = act.get_MT_rmax() * bonds[i].nominalLengthInAngstroms;
            bonds[i].hill_co = act.get_MT_hill_co();
            
            
            
            bonds[i+act_MT_bonds].type = 7;
            bonds[i+act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+act_MT_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+act_MT_bonds].nominalLengthInAngstroms=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            bonds[i+act_MT_bonds].stiffnessInKcalPerAngstrom2 = act.get_MT_k1();
            bonds[i+act_MT_bonds].r_min = act.get_MT_rmin() * bonds[i+act_MT_bonds].nominalLengthInAngstroms;
            bonds[i+act_MT_bonds].r_max = bonds[i+act_MT_bonds].nominalLengthInAngstroms;
            
            
            
            bonds[i+2*act_MT_bonds].type = 7;
            bonds[i+2*act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+2*act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+2*act_MT_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+2*act_MT_bonds].nominalLengthInAngstroms=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            bonds[i+2*act_MT_bonds].stiffnessInKcalPerAngstrom2 = act.get_MT_k2();
            bonds[i+2*act_MT_bonds].r_min = bonds[i+2*act_MT_bonds].nominalLengthInAngstroms;
            bonds[i+2*act_MT_bonds].r_max = act.get_MT_rmax() * bonds[i+2*act_MT_bonds].nominalLengthInAngstroms;
            
            
            
            
            
            bonds[i+3*act_MT_bonds].type = act.get_MT_spring_model();
            bonds[i+3*act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+3*act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+3*act_MT_bonds].class_label = act.get_label() + act.get_label();
            switch (bonds[i+3*act_MT_bonds].type) {
                    //FENE
                case 1:
                    //                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                    //                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                    break;
                    //Harmonic
                case 2:
                    bonds[i+3*act_MT_bonds].nominalLengthInAngstroms=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                    bonds[i+3*act_MT_bonds].stiffnessInKcalPerAngstrom2=act.get_MT_spring_stiffness_coefficient();
                    break;
                    
                    //Kelvin-Voigt
                case 4:
                    bonds[i+3*act_MT_bonds].nominalLengthInAngstroms=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                    bonds[i+3*act_MT_bonds].stiffnessInKcalPerAngstrom2=act.get_MT_spring_stiffness_coefficient();
                    bonds[i+3*act_MT_bonds].dampInKcalPsPerAngstrom2=act.get_kelvin_damping_coefficient();
                    break;
                    
            }
            
            
        }
    
    
    
    if(bonds[0].type==2){
        cout<<"Actin bond potential: Harmonic "<<endl;
    }
    
    if(bonds[0].type == 4){
        cout<<"Actin bond potential: Kelvin-Voigt "<<endl;
        //cout<<"spring coeficient (KJ per Nanometer4) ="<< act.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
        //cout<<"bending coeficient (KJ per Nanometer2)="<<mem.get_bending_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    
    if(bonds[act_num_bonds].type == 6){
        cout<<"Actin Contractile element with constant force "<< act.get_contractile_force()<< endl;
    }
    
    if(bonds[2*act_num_bonds].type == 7){
        cout<<"Actin Contractile element with k1 and k2 "<< act.get_contractile_k1()<< " and " << act.get_contractile_k2()<< endl;
        cout<<"Actin Contractile element with rmin and rmax "<< act.get_contractile_rmin() * act.get_avg_node_dist()<< " and " << act.get_contractile_rmax() * act.get_avg_node_dist()<< endl;
    }
    
    cout<<endl;
    
    return bonds;
}



Bonds* convert_ActMem_bond_info_to_openmm(Actin act, int k) {
    const int act_num_bonds = act.return_num_of_actin_membrane_shared_nodes(k);
    Bonds* bonds = new Bonds[act_num_bonds];
    double stiffness = 1000000;
    double nominal_length = 0;
    int type = 2;
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<act_num_bonds; i++) {
        bonds[i].type = type;
        bonds[i].atoms[0]=act.return_ActMem_shared_act_atom(k, i);
        bonds[i].atoms[1]=act.return_ActMem_shared_mem_atom(k, i);
        bonds[i].nominalLengthInAngstroms=nominal_length;
        bonds[i].stiffnessInKcalPerAngstrom2=stiffness;
        bonds[i].class_label = act.get_label() + GenConst::Membrane_label+std::to_string(k);
    }
    //std::cout << "convert" << bonds[1].class_label << '\n';
    return bonds;
}


Bonds* convert_ECM_bond_info_to_openmm(ECM ecm , MyAtomInfo* atoms) {
    const int ecm_num_bonds = ecm.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[ecm_num_bonds];
//    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
//    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<ecm_num_bonds; i++) {
        bonds[i].type = ecm.get_spring_model();
        bonds[i].atoms[0]=ecm.get_node_pair(i, 0);
        bonds[i].atoms[1]=ecm.get_node_pair(i, 1);
        bonds[i].class_label = ecm.get_label() + ecm.get_label();
        
        
        double dist;
        dist = sqrt((atoms[ bonds[i].atoms[0] ].initPosInAng[0] - atoms[ bonds[i].atoms[1] ].initPosInAng[0]) * (atoms[ bonds[i].atoms[0] ].initPosInAng[0] - atoms[ bonds[i].atoms[1] ].initPosInAng[0]) + (atoms[ bonds[i].atoms[0] ].initPosInAng[1] - atoms[ bonds[i].atoms[1] ].initPosInAng[1]) * (atoms[ bonds[i].atoms[0] ].initPosInAng[1] - atoms[ bonds[i].atoms[1] ].initPosInAng[1]) + (atoms[ bonds[i].atoms[0] ].initPosInAng[2] - atoms[ bonds[i].atoms[1] ].initPosInAng[2]) * (atoms[ bonds[i].atoms[0] ].initPosInAng[2] - atoms[ bonds[i].atoms[1] ].initPosInAng[2])) ;
        
        
        switch (bonds[i].type) {
                //FENE
            case 1:
                ecm.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                bonds[i].stiffnessInKcalPerAngstrom2=ecm.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                //bonds[i].nominalLengthInAngstroms=ecm.get_avg_node_dist();
                bonds[i].nominalLengthInAngstroms=dist;
                
                
                bonds[i].stiffnessInKcalPerAngstrom2=ecm.get_spring_stiffness_coefficient() +0.5*(atoms[bonds[i].atoms[0]].initPosInAng[0] + atoms[bonds[i].atoms[1]].initPosInAng[0]) * ecm.get_stiffness_gradient_x()
                +0.5*(atoms[bonds[i].atoms[0]].initPosInAng[1] + atoms[bonds[i].atoms[1]].initPosInAng[1]) * ecm.get_stiffness_gradient_y()
                +0.5*(atoms[bonds[i].atoms[0]].initPosInAng[2] + atoms[bonds[i].atoms[1]].initPosInAng[2]) * ecm.get_stiffness_gradient_z();
               
                break;
                
                
        }
        
        
    }

    return bonds;
}


Bonds* convert_Chromatin_bond_info_to_openmm(Chromatin chromo) {
    const int chromo_num_bonds = chromo.get_num_of_nodes()-1;
    Bonds* bonds = new Bonds[chromo_num_bonds];
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<chromo_num_bonds; i++) {
        bonds[i].type = chromo.get_spring_model();
        bonds[i].atoms[0]=i;
        bonds[i].atoms[1]=i+1;
        bonds[i].class_label = chromo.get_label() + chromo.get_label();
        switch (bonds[i].type) {
                //FENE
            case 1:
                //                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                //                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                break;
                //Harmonic
            case 2:
                bonds[i].nominalLengthInAngstroms=2*chromo.get_node_radius();
                bonds[i].stiffnessInKcalPerAngstrom2=chromo.get_spring_stiffness_coefficient();
                break;
                
                
        }
        
        
    }
    
    return bonds;
}
