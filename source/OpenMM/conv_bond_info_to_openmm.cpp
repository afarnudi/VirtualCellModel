#include "OpenMM_funcs.hpp"

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = mem.get_spring_model();
        bonds[i].atoms[0]=mem.get_node_pair(i, 0);
        bonds[i].atoms[1]=mem.get_node_pair(i, 1);
        bonds[i].class_label = mem.get_label() + mem.get_label();
//        if (mem.get_new_node_radius()!=-1) {
//            bonds[i].nominalLengthInNm=mem.get_node_pair_distance(i);
//        } else {
//            bonds[i].nominalLengthInNm=2*mem.get_new_node_radius();
//        }
        bonds[i].nominalLengthInNm=mem.get_node_pair_distance(i);
        
        switch (bonds[i].type) {
                //FENE
            case 1:
                mem.set_FENE_param_2(bonds[i].FENE_lmininNm, bonds[i].FENE_lmaxinNm, bonds[i].FENE_epsilon, bonds[i].FENE_k);
                break;
                //Harmonic
            case 2:
                bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
                break;
                //HarmonicX4
            case 3:
                bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
                break;
                //Kelvin-Voigt
            case 4:
                bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
                break;
                
            case 5:
                bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
                break;
        }
        
        
    }
    if(bonds[0].type==2){
        cout<<" Harmonic "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<mem.get_spring_stiffness_coefficient() <<endl;
    }
    
    if (bonds[0].type == 1) {
        cout<< " FENE"<<endl;
//        cout<<"attraction coeficient (KJ . Nm^-2 . mol^-1 ) = "<<mem.get_spring_stiffness_coefficient() <<endl;
//        cout<<"bending coeficient (KJ . mol^-1 ) = "<<mem.get_bending_stiffness_coefficient() <<endl;
//        cout<<"lmin(Nm)\tlmax(Nm)\tK_FENE   \tle1(Nm)\n"<<bonds[0].FENE_lmininNm<<"\t"<<bonds[0].FENE_lmaxinNm<<"\t"<<bonds[0].K_FENE<<"\t"<<bonds[0].FENE_le1inNm<<endl;
    }
    
    if(bonds[0].type == 3){
        cout<<" X4Harmonic "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-4 . mol^-1 ) = "<<mem.get_spring_stiffness_coefficient()  <<endl;
    }
    
    if(bonds[0].type == 4){
        cout<<" Kelvin-Voigt "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-2 . mol^-1 ) = " <<mem.get_spring_stiffness_coefficient() <<endl;
    }
    if(bonds[0].type==5){
        cout<<" realHarmonic "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-2 . mol^-1 ) = "<<mem.get_spring_stiffness_coefficient() <<endl;
    }
    
    
    
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
                // bonds[i].nominalLengthInNm=act.get_avg_node_dist();
                bonds[i].nominalLengthInNm=act.get_act_relaxlength(i) * act.get_act_r0factor();
                bonds[i].stiffnessInKJPerNm2=act.get_spring_stiffness_coefficient();

                break;
                
                //Kelvin-Voigt
            case 4:
                //bonds[i].nominalLengthInNm=act.get_avg_node_dist();
                bonds[i].nominalLengthInNm=act.get_act_relaxlength(i) * act.get_act_r0factor();
                bonds[i].stiffnessInKJPerNm2=act.get_spring_stiffness_coefficient();
                bonds[i].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
                
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
        
        //bonds[i].nominalLengthInNm=act.get_avg_node_dist();
        bonds[i].nominalLengthInNm=act.get_act_relaxlength(i-act_num_bonds) * act.get_act_r0factor();
        bonds[i].F0 = act.get_contractile_force();
        bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInNm ;
        bonds[i].r_max = act.get_contractile_rmax() * bonds[i].nominalLengthInNm;
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
        
        //bonds[i].nominalLengthInNm=act.get_avg_node_dist();
        bonds[i].nominalLengthInNm=act.get_act_relaxlength(i-2*act_num_bonds) * act.get_act_r0factor();
        bonds[i].stiffnessInKJPerNm2 = act.get_contractile_k1();
        bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInNm;
        bonds[i].r_max = bonds[i].nominalLengthInNm;
        //bonds[i].r_max = act.get_avg_node_dist();
        
        bonds[i+act_num_bonds].nominalLengthInNm=act.get_act_relaxlength(i-2*act_num_bonds) * act.get_act_r0factor();
        bonds[i+act_num_bonds].stiffnessInKJPerNm2 = act.get_contractile_k2();
        bonds[i+act_num_bonds].r_min = bonds[i+act_num_bonds].nominalLengthInNm;
        bonds[i+act_num_bonds].r_max = act.get_contractile_rmax() * bonds[i+act_num_bonds].nominalLengthInNm;
    }
    
    
//    int n1=0;
//    int n2=0;
    
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
        
        bonds[i].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        
        
        
        bonds[i].F0 = act.get_abp_force();
        
//        if( (0.5*(atoms[bonds[i].atoms[0]].initPosInNm[0] + atoms[bonds[i].atoms[1]].initPosInNm[0]) ) >0  )
//        {
//            bonds[i].F0 =  act.get_abp_force();
//            n1++;
//        }
//        else
//        {
//            bonds[i].F0 = 10*act.get_abp_force();
//            n2++;
//        }
        
        
        bonds[i].r_min = act.get_abp_rmin() * bonds[i].nominalLengthInNm;
        bonds[i].r_max = act.get_abp_rmax() * bonds[i].nominalLengthInNm;
        bonds[i].hill_co = act.get_abp_hill_co();
        
        
        
        bonds[i+act_abp_bonds].type = 7;
        bonds[i+act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i+act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i+act_abp_bonds].class_label = act.get_label() + act.get_label();
        
        bonds[i+act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        bonds[i+act_abp_bonds].stiffnessInKJPerNm2 = act.get_abp_k1();
        bonds[i+act_abp_bonds].r_min = act.get_abp_rmin() * bonds[i+act_abp_bonds].nominalLengthInNm;
        bonds[i+act_abp_bonds].r_max = bonds[i+act_abp_bonds].nominalLengthInNm;
        
        
        
        bonds[i+2*act_abp_bonds].type = 7;
        bonds[i+2*act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
        bonds[i+2*act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
        bonds[i+2*act_abp_bonds].class_label = act.get_label() + act.get_label();
        
        bonds[i+2*act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
        bonds[i+2*act_abp_bonds].stiffnessInKJPerNm2 = act.get_abp_k2();
        bonds[i+2*act_abp_bonds].r_min = bonds[i+2*act_abp_bonds].nominalLengthInNm;
        bonds[i+2*act_abp_bonds].r_max = act.get_abp_rmax() * bonds[i+2*act_abp_bonds].nominalLengthInNm;
        
        
        
        
        
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
                bonds[i+3*act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKJPerNm2=act.get_abp_spring_stiffness_coefficient();
                break;
                
                //Kelvin-Voigt
            case 4:
                bonds[i+3*act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKJPerNm2=act.get_abp_spring_stiffness_coefficient();
                bonds[i+3*act_abp_bonds].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
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
            
            bonds[i].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            
            
            
            bonds[i].F0 = act.get_MT_force();
            
//            if( (0.5*(atoms[bonds[i].atoms[0]].initPosInNm[0] + atoms[bonds[i].atoms[1]].initPosInNm[0]) ) >0  )
//            {
//                bonds[i].F0 =  act.get_MT_force();
//                n1++;
//            }
//            else
//            {
//                bonds[i].F0 = 5*act.get_MT_force();
//                n2++;
//            }
            
            
            bonds[i].r_min = act.get_MT_rmin() * bonds[i].nominalLengthInNm;
            bonds[i].r_max = act.get_MT_rmax() * bonds[i].nominalLengthInNm;
            bonds[i].hill_co = act.get_MT_hill_co();
            
            
            
            bonds[i+act_MT_bonds].type = 7;
            bonds[i+act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+act_MT_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            bonds[i+act_MT_bonds].stiffnessInKJPerNm2 = act.get_MT_k1();
            bonds[i+act_MT_bonds].r_min = act.get_MT_rmin() * bonds[i+act_MT_bonds].nominalLengthInNm;
            bonds[i+act_MT_bonds].r_max = bonds[i+act_MT_bonds].nominalLengthInNm;
            
            
            
            bonds[i+2*act_MT_bonds].type = 7;
            bonds[i+2*act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+2*act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+2*act_MT_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+2*act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            bonds[i+2*act_MT_bonds].stiffnessInKJPerNm2 = act.get_MT_k2();
            bonds[i+2*act_MT_bonds].r_min = bonds[i+2*act_MT_bonds].nominalLengthInNm;
            bonds[i+2*act_MT_bonds].r_max = act.get_MT_rmax() * bonds[i+2*act_MT_bonds].nominalLengthInNm;
            
            
            
            
            
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
                    bonds[i+3*act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                    bonds[i+3*act_MT_bonds].stiffnessInKJPerNm2=act.get_MT_spring_stiffness_coefficient();
                    break;
                    
                    //Kelvin-Voigt
                case 4:
                    bonds[i+3*act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                    bonds[i+3*act_MT_bonds].stiffnessInKJPerNm2=act.get_MT_spring_stiffness_coefficient();
                    bonds[i+3*act_MT_bonds].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
                    break;
                    
            }
            
            
        }
    
    
    
    if(bonds[0].type==2){
        cout<<" Harmonic "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<act.get_spring_stiffness_coefficient() <<endl;
    }
    
    if(bonds[0].type == 4){
        cout<<" Kelvin-Voigt "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-4) = "<< act.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    
    if(bonds[act_num_bonds].type == 6){
        cout<<"Contractile element with constant force "<<endl;
        cout<<"\tContractile force? (?units?) = "<<act.get_contractile_force()<< endl;
    }
    
    if(bonds[2*act_num_bonds].type == 7){
        cout<<" Contractile elements"<<endl;
        cout<<"\tk1 and k2 ? (?units?) = "<< act.get_contractile_k1()<< " , " << act.get_contractile_k2()<< endl;
        cout<<"\trmin and rmax (Nm) = "<< act.get_contractile_rmin() * act.get_avg_node_dist()<< " , " << act.get_contractile_rmax() * act.get_avg_node_dist()<< endl;
    }
    
    
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
        bonds[i].nominalLengthInNm=nominal_length;
        bonds[i].stiffnessInKJPerNm2=stiffness;
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
        dist = sqrt((atoms[ bonds[i].atoms[0] ].initPosInNm[0] - atoms[ bonds[i].atoms[1] ].initPosInNm[0]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[0] - atoms[ bonds[i].atoms[1] ].initPosInNm[0]) + (atoms[ bonds[i].atoms[0] ].initPosInNm[1] - atoms[ bonds[i].atoms[1] ].initPosInNm[1]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[1] - atoms[ bonds[i].atoms[1] ].initPosInNm[1]) + (atoms[ bonds[i].atoms[0] ].initPosInNm[2] - atoms[ bonds[i].atoms[1] ].initPosInNm[2]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[2] - atoms[ bonds[i].atoms[1] ].initPosInNm[2])) ;
        
        
        switch (bonds[i].type) {
                //FENE
            case 1:
//                ecm.set_FENE_param(bonds[i].FENE_lmininNm, bonds[i].FENE_lmaxinNm, bonds[i].FENE_epsilon, bonds[i].FENE_k);
//                bonds[i].stiffnessInKJPerNm2=ecm.get_spring_stiffness_coefficient();
                cout<<TWWARN<<"FENE not set"<<TRESET<<endl;
                exit(0);
                break;
                //Harmonic
            case 2:
                //bonds[i].nominalLengthInNm=ecm.get_avg_node_dist();
                bonds[i].nominalLengthInNm=dist;
                
                
                bonds[i].stiffnessInKJPerNm2=ecm.get_spring_stiffness_coefficient() +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[0] + atoms[bonds[i].atoms[1]].initPosInNm[0]) * ecm.get_stiffness_gradient_x()
                +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[1] + atoms[bonds[i].atoms[1]].initPosInNm[1]) * ecm.get_stiffness_gradient_y()
                +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[2] + atoms[bonds[i].atoms[1]].initPosInNm[2]) * ecm.get_stiffness_gradient_z();
               
                break;
                
                
        }
        if(bonds[0].type==2){
            cout<<" Harmonic "<<endl;
            cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<ecm.get_spring_stiffness_coefficient() <<endl;
        }
        
    }

    return bonds;
}


Bonds* convert_Chromatin_bond_info_to_openmm(Chromatin chromo) {
    const int chromo_num_bonds = chromo.get_num_of_bonds();
    Bonds* bonds = new Bonds[chromo_num_bonds];
    
    
    cout<<"chromo_num_bonds = "<<chromo_num_bonds<<endl;
    if (chromo_num_bonds != 0) {
        int num_of_real_bonds = chromo.get_num_of_real_site()-1;
            int num_virtual_sites_per_bond= chromo.get_num_of_virtual_sites_per_bond();
            //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
            //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
            for (int i=0; i<num_of_real_bonds; i++) {
                bonds[i].type = chromo.get_spring_model();
                bonds[i].atoms[0]=i*(num_virtual_sites_per_bond+1);
                bonds[i].atoms[1]=(i+1)*(num_virtual_sites_per_bond+1);
        //        cout<<"atom bon 1 = "<<i*(num_virtual_sites_per_bond+1)<<endl;
        //        cout<<"atom bon 2 = "<<(i+1)*(num_virtual_sites_per_bond+1)<<endl;
                bonds[i].class_label = chromo.get_label() + chromo.get_label();
                switch (bonds[i].type) {
                        //FENE
                    case 1:
                        //                act.set_FENE_param(bonds[i].FENE_le0, bonds[i].FENE_le1, bonds[i].FENE_lmin, bonds[i].FENE_lmax);
                        //                bonds[i].stiffnessInKcalPerAngstrom2=act.get_spring_stiffness_coefficient();
                        break;
                        //Harmonic
                    case 2:
                        bonds[i].nominalLengthInNm=chromo.get_bond_length();
                        bonds[i].stiffnessInKJPerNm2=chromo.get_spring_stiffness_coefficient();
                        break;
                        
                        
                }
                
                
            }
        if(bonds[0].type==2){
            cout<<" Harmonic "<<endl;
            cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<chromo.get_spring_stiffness_coefficient() <<endl;
        }
        vector<vector<int> > virtual_bond_list = chromo.get_virtual_bonds();
        
        int list_counter=0;
        for (int i=num_of_real_bonds+1; i<chromo_num_bonds; i++) {
            
            bonds[i].type = -2;
            bonds[i].atoms[0]=virtual_bond_list[list_counter][0];
            bonds[i].atoms[1]=virtual_bond_list[list_counter][1];
            bonds[i].class_label = chromo.get_label() + chromo.get_label();
            bonds[i].nominalLengthInNm=0;
            bonds[i].stiffnessInKJPerNm2=0;
        }
    }
    return bonds;
}
