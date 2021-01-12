#include "OpenMM_funcs.hpp"

Bonds* convert_membrane_bond_info_to_openmm(Membrane mem) {
    const int mem_num_bonds = mem.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[mem_num_bonds];
    
    bool fenepotential =false;
    bool harmonicpotential = false;
    bool gompperpotential = false;
    
    
    for (int i=0; i<mem_num_bonds; i++) {
        bonds[i].type = mem.get_spring_model();
        bonds[i].atoms[0]=mem.get_node_pair(i, 0);
        bonds[i].atoms[1]=mem.get_node_pair(i, 1);
        bonds[i].class_label = mem.get_label() + mem.get_label();
        bonds[i].nominalLengthInNm=mem.get_node_pair_Nominal_Length_in_Nm(i);
        
        if (bonds[i].type == potentialModelIndex.Model["FENE"]) {
            fenepotential =true;
            bonds[i].FENER0inNm = 1.5*bonds[i].nominalLengthInNm;
            bonds[i].k_FENE_inKJpermol = 30*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature/(bonds[i].nominalLengthInNm*bonds[i].nominalLengthInNm);
            bonds[i].epsilon_FENE_inKJpermol = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
        } else if (bonds[i].type == potentialModelIndex.Model["Harmonic"]){
            harmonicpotential=true;
            bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["HarmonicX4"]){
            bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["Kelvin-Voigt"]){
            bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["RealHarmonic"]){
            bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["Gompper"]){
            harmonicpotential=true;
            vector<double> gompperparamslminlc1lc0lmax = mem.get_gompper_params_lminlc1lc0lmax();
            bonds[i].stiffnessInKJPerNm2=mem.get_spring_stiffness_coefficient();
            bonds[i].gompperlmin=gompperparamslminlc1lc0lmax[0];
            bonds[i].gompperlc1=gompperparamslminlc1lc0lmax[1];
            bonds[i].gompperlc0=gompperparamslminlc1lc0lmax[2];
            bonds[i].gompperlmax=gompperparamslminlc1lc0lmax[3];
        }
        
    }
    if(harmonicpotential){
        cout<<" Harmonic "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<mem.get_spring_stiffness_coefficient() <<endl;
    }
    
    if (fenepotential) {
        cout<< " FENE"<<endl;
    }
    
    if(bonds[0].type == potentialModelIndex.Model["HarmonicX4"]){
        cout<<" X4Harmonic "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-4 . mol^-1 ) = "<<mem.get_spring_stiffness_coefficient()  <<endl;
    }
    
    if(bonds[0].type == potentialModelIndex.Model["Kelvin-Voigt"]){
        cout<<" Kelvin-Voigt "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-2 . mol^-1 ) = " <<mem.get_spring_stiffness_coefficient() <<endl;
    }
    if(bonds[0].type == potentialModelIndex.Model["RealHarmonic"]){
        cout<<" realHarmonic "<<endl;
        cout<<"\tCoeficient (KJ . Nm^-2 . mol^-1 ) = "<<mem.get_spring_stiffness_coefficient() <<endl;
    }
    
    if(gompperpotential){
        cout<<" Gompper "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<mem.get_spring_stiffness_coefficient() <<endl;
        cout<<"\lmin,   lc1,   lc0,   lmax (Nm): \n\t" <<mem.get_spring_stiffness_coefficient() <<endl;
    }
    
    return bonds;
}

Bonds* convert_Actin_bond_info_to_openmm(Actin act,MyAtomInfo* atoms) {
    const int act_num_bonds = act.get_num_of_node_pairs();
    const int act_abp_bonds = act.get_num_of_abp_pairs();
    const int act_MT_bonds = act.get_num_of_MT_pairs();
    
    bool harmonicpotential= false;
    bool KelvinVoigtpotential = false;
    bool Contractilepotential = false;
    bool Harmonic_minmaxpotential = false;
    bool fenepotential=false;
    
    Bonds* bonds = new Bonds[act_num_bonds];
//    Bonds* bonds = new Bonds[4*act_num_bonds + 4*act_abp_bonds + 4*act_MT_bonds];
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<act_num_bonds; i++) {
        bonds[i].type = act.get_spring_model();
        bonds[i].atoms[0]=act.get_node_pair(i, 0);
        bonds[i].atoms[1]=act.get_node_pair(i, 1);
        bonds[i].class_label = act.get_label() + act.get_label();
        
        if (bonds[i].type == potentialModelIndex.Model["FENE"]) {
            fenepotential =true;
            bonds[i].FENER0inNm = 1.5*bonds[i].nominalLengthInNm;
            bonds[i].k_FENE_inKJpermol = 30*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature/(bonds[i].nominalLengthInNm*bonds[i].nominalLengthInNm);
            bonds[i].epsilon_FENE_inKJpermol = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
        } else if (bonds[i].type == potentialModelIndex.Model["Harmonic"]){
            harmonicpotential= true;
            bonds[i].nominalLengthInNm=act.get_node_pair_Nominal_Length_in_Nm(i) * act.get_act_r0factor();
            bonds[i].stiffnessInKJPerNm2=act.get_spring_stiffness_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["Kelvin-Voigt"]){
            KelvinVoigtpotential = true;
            bonds[i].nominalLengthInNm=act.get_node_pair_Nominal_Length_in_Nm(i) * act.get_act_r0factor();
            bonds[i].stiffnessInKJPerNm2=act.get_spring_stiffness_coefficient();
            bonds[i].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
        } else if (bonds[i].type == potentialModelIndex.Model["Contractile"]){
            Contractilepotential = true;
        }
        
        
        
    }
    
    
    //Contractile element parallel to bond spring
    if (Contractilepotential) {
        for (int i=act_num_bonds; i<(2*act_num_bonds); i++) {
            int model = act.get_contractile_model();
            switch (model) {
                case 1:
                    bonds[i].type = potentialModelIndex.Model["Contractile"];
                    Contractilepotential = true;
                    break;
                case 2:
                    bonds[i].type = potentialModelIndex.Model["hill"];
                    break;
                case 3:
                    bonds[i].type = potentialModelIndex.Model["KFs"];
                    break;
            }
            bonds[i].k_F0 = 4;
            //bonds[i].type = 6;
            bonds[i].atoms[0]=act.get_node_pair(i-act_num_bonds, 0);
            bonds[i].atoms[1]=act.get_node_pair(i-act_num_bonds, 1);
            bonds[i].class_label = act.get_label() + act.get_label();
            
            //bonds[i].nominalLengthInNm=act.get_avg_node_dist();
            bonds[i].nominalLengthInNm=act.get_node_pair_Nominal_Length_in_Nm(i-act_num_bonds) * act.get_act_r0factor();
            bonds[i].F0 = act.get_contractile_force();
            bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInNm ;
            bonds[i].r_max = act.get_contractile_rmax() * bonds[i].nominalLengthInNm;
            //bonds[i].hill_co = act.get_hill_co();
            bonds[i].hill_co = act.get_contractile_hill_co();
        }
        
        //Contractile spring
        for (int i=2*act_num_bonds; i<(3*act_num_bonds); i++) {
            bonds[i].type = potentialModelIndex.Model["Harmonic"];
            bonds[i+act_num_bonds].type = potentialModelIndex.Model["Harmonic_minmax"];
            bonds[i].atoms[0]=act.get_node_pair(i-2*act_num_bonds, 0);
            bonds[i].atoms[1]=act.get_node_pair(i-2*act_num_bonds, 1);
            bonds[i+act_num_bonds].atoms[0]=act.get_node_pair(i-2*act_num_bonds, 0);
            bonds[i+act_num_bonds].atoms[1]=act.get_node_pair(i-2*act_num_bonds, 1);
            bonds[i].class_label = act.get_label() + act.get_label();
            bonds[i+act_num_bonds].class_label = act.get_label() + act.get_label();
            
            //bonds[i].nominalLengthInNm=act.get_avg_node_dist();
            bonds[i].nominalLengthInNm=act.get_node_pair_Nominal_Length_in_Nm(i-2*act_num_bonds) * act.get_act_r0factor();
            bonds[i].stiffnessInKJPerNm2 = act.get_contractile_k1();
            bonds[i].r_min = act.get_contractile_rmin() * bonds[i].nominalLengthInNm;
            bonds[i].r_max = bonds[i].nominalLengthInNm;
            //bonds[i].r_max = act.get_avg_node_dist();
            
            bonds[i+act_num_bonds].nominalLengthInNm=act.get_node_pair_Nominal_Length_in_Nm(i-2*act_num_bonds) * act.get_act_r0factor();
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
                    bonds[i].type = potentialModelIndex.Model["Contractile"];
                    Contractilepotential = true;
                    break;
                case 2:
                    bonds[i].type = potentialModelIndex.Model["hill"];
                    break;
                case 3:
                    bonds[i].type = potentialModelIndex.Model["KFs"];
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
            
            
            
            bonds[i+act_abp_bonds].type = potentialModelIndex.Model["Harmonic_minmax"];
            bonds[i+act_abp_bonds].atoms[0]=act.get_abp_pair(i-4*act_num_bonds, 0);
            bonds[i+act_abp_bonds].atoms[1]=act.get_abp_pair(i-4*act_num_bonds, 1);
            bonds[i+act_abp_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
            bonds[i+act_abp_bonds].stiffnessInKJPerNm2 = act.get_abp_k1();
            bonds[i+act_abp_bonds].r_min = act.get_abp_rmin() * bonds[i+act_abp_bonds].nominalLengthInNm;
            bonds[i+act_abp_bonds].r_max = bonds[i+act_abp_bonds].nominalLengthInNm;
            
            
            
            bonds[i+2*act_abp_bonds].type = potentialModelIndex.Model["Harmonic_minmax"];
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
            
            if (bonds[i+3*act_abp_bonds].type == potentialModelIndex.Model["Harmonic"]) {
                bonds[i+3*act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKJPerNm2=act.get_abp_spring_stiffness_coefficient();
            } else if (bonds[i+3*act_abp_bonds].type == potentialModelIndex.Model["Kelvin-Voigt"]){
                bonds[i+3*act_abp_bonds].nominalLengthInNm=act.get_abp_r0factor() * act.get_abp_relaxlength(i-4*act_num_bonds) ;
                bonds[i+3*act_abp_bonds].stiffnessInKJPerNm2=act.get_abp_spring_stiffness_coefficient();
                bonds[i+3*act_abp_bonds].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
            }
        }
        
        
        
        
        
        
        
        for (int i=4*act_num_bonds+ 4*act_abp_bonds; i<(4*act_num_bonds + 4*act_abp_bonds + act_MT_bonds); i++) {
            
            int model = act.get_MT_model();
            switch (model) {
                case 1:
                    bonds[i].type = potentialModelIndex.Model["Contractile"];
                    Contractilepotential = true;
                    break;
                case 2:
                    bonds[i].type = potentialModelIndex.Model["hill"];
                    break;
                case 3:
                    bonds[i].type = potentialModelIndex.Model["KFs"];
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
            
            
            
            bonds[i+act_MT_bonds].type = potentialModelIndex.Model["Harmonic_minmax"];
            bonds[i+act_MT_bonds].atoms[0]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 0);
            bonds[i+act_MT_bonds].atoms[1]=act.get_MT_pair(i-4*act_num_bonds - 4*act_abp_bonds, 1);
            bonds[i+act_MT_bonds].class_label = act.get_label() + act.get_label();
            
            bonds[i+act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
            bonds[i+act_MT_bonds].stiffnessInKJPerNm2 = act.get_MT_k1();
            bonds[i+act_MT_bonds].r_min = act.get_MT_rmin() * bonds[i+act_MT_bonds].nominalLengthInNm;
            bonds[i+act_MT_bonds].r_max = bonds[i+act_MT_bonds].nominalLengthInNm;
            
            
            
            bonds[i+2*act_MT_bonds].type = potentialModelIndex.Model["Harmonic_minmax"];
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
            
            if (bonds[i+3*act_MT_bonds].type == potentialModelIndex.Model["Harmonic"]) {
                bonds[i+3*act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                bonds[i+3*act_MT_bonds].stiffnessInKJPerNm2=act.get_MT_spring_stiffness_coefficient();
            } else if (bonds[i+3*act_MT_bonds].type == potentialModelIndex.Model["Kelvin-Voigt"]){
                bonds[i+3*act_MT_bonds].nominalLengthInNm=act.get_MT_r0factor() * act.get_MT_relaxlength(i-4*act_num_bonds - 4*act_abp_bonds) ;
                bonds[i+3*act_MT_bonds].stiffnessInKJPerNm2=act.get_MT_spring_stiffness_coefficient();
                bonds[i+3*act_MT_bonds].dampInKJPsPerNm2=act.get_kelvin_damping_coefficient();
            }
            
            
            
            
        }
        
        
    }
    
    
    if(harmonicpotential){
        cout<<" Harmonic "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<act.get_spring_stiffness_coefficient() <<endl;
    }
    
    if(KelvinVoigtpotential){
        cout<<" Kelvin-Voigt "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-4) = "<< act.get_spring_stiffness_coefficient() * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm<<endl;
    }
    
    if(Contractilepotential){
        cout<<"Contractile element with constant force "<<endl;
        cout<<"\tContractile force? (?units?) = "<<act.get_contractile_force()<< endl;
    }
    
    if(bonds[2*act_num_bonds].type == potentialModelIndex.Model["Harmonic_minmax"]){
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
    int type = potentialModelIndex.Model["Harmonic"];
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<act_num_bonds; i++) {
        bonds[i].type = type;
        bonds[i].atoms[0]=act.return_ActMem_shared_act_atom(k, i);
        bonds[i].atoms[1]=act.return_ActMem_shared_mem_atom(k, i);
        bonds[i].nominalLengthInNm=nominal_length;
        bonds[i].stiffnessInKJPerNm2=stiffness;
        bonds[i].class_label = act.get_label() + generalParameters.Membrane_label+std::to_string(k);
    }
    //std::cout << "convert" << bonds[1].class_label << '\n';
    return bonds;
}


Bonds* convert_ECM_bond_info_to_openmm(ECM ecm , MyAtomInfo* atoms) {
    const int ecm_num_bonds = ecm.get_num_of_node_pairs();
    Bonds* bonds = new Bonds[ecm_num_bonds];
    bool harmonicpotential= false;
    //    cout<<"ecm.get_spring_model()  "<<ecm.get_spring_model()<<endl;
    //    cout<<"ecm.get_num_of_node_pairs()()  "<<ecm.get_num_of_node_pairs()<<endl;
    for (int i=0; i<ecm_num_bonds; i++) {
        bonds[i].type = ecm.get_spring_model();
        bonds[i].atoms[0]=ecm.get_node_pair(i, 0);
        bonds[i].atoms[1]=ecm.get_node_pair(i, 1);
        bonds[i].class_label = ecm.get_label() + ecm.get_label();
        
        
        double dist;
        dist = sqrt((atoms[ bonds[i].atoms[0] ].initPosInNm[0] - atoms[ bonds[i].atoms[1] ].initPosInNm[0]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[0] - atoms[ bonds[i].atoms[1] ].initPosInNm[0]) + (atoms[ bonds[i].atoms[0] ].initPosInNm[1] - atoms[ bonds[i].atoms[1] ].initPosInNm[1]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[1] - atoms[ bonds[i].atoms[1] ].initPosInNm[1]) + (atoms[ bonds[i].atoms[0] ].initPosInNm[2] - atoms[ bonds[i].atoms[1] ].initPosInNm[2]) * (atoms[ bonds[i].atoms[0] ].initPosInNm[2] - atoms[ bonds[i].atoms[1] ].initPosInNm[2])) ;
        
        if (bonds[i].type == potentialModelIndex.Model["FENE"]) {
            cout<<TWWARN<<"FENE not set"<<TRESET<<endl;
            exit(0);
        } else if (bonds[i].type == potentialModelIndex.Model["Harmonic"]){
            harmonicpotential=true;
            bonds[i].nominalLengthInNm=dist;
            
            
            bonds[i].stiffnessInKJPerNm2=ecm.get_spring_stiffness_coefficient() +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[0] + atoms[bonds[i].atoms[1]].initPosInNm[0]) * ecm.get_stiffness_gradient_x()
            +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[1] + atoms[bonds[i].atoms[1]].initPosInNm[1]) * ecm.get_stiffness_gradient_y()
            +0.5*(atoms[bonds[i].atoms[0]].initPosInNm[2] + atoms[bonds[i].atoms[1]].initPosInNm[2]) * ecm.get_stiffness_gradient_z();
        }
        
    }
    
    if(harmonicpotential){
        cout<<" Harmonic "<<endl;
        cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<ecm.get_spring_stiffness_coefficient() <<endl;
    }
    
    return bonds;
}


Bonds* convert_Chromatin_bond_info_to_openmm(Chromatin chromo) {
    const int chromo_num_bonds = chromo.get_num_of_bonds();
    Bonds* bonds = new Bonds[chromo_num_bonds];
    
    bool harmonicpotential =false;
    bool fenepotential =false;
    //    cout<<"\nchromo_num_bonds = "<<chromo_num_bonds<<endl;
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
            
            if (bonds[i].type == potentialModelIndex.Model["FENE"]) {
                fenepotential =true;
                bonds[i].nominalLengthInNm=chromo.get_bond_nominal_length(i);
                bonds[i].FENER0inNm = 1.5*bonds[i].nominalLengthInNm;
                bonds[i].k_FENE_inKJpermol = 30*generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature/(bonds[i].nominalLengthInNm*bonds[i].nominalLengthInNm);
                bonds[i].epsilon_FENE_inKJpermol = generalParameters.BoltzmannKJpermolkelvin*generalParameters.temperature;
                
            } else if (bonds[i].type == potentialModelIndex.Model["Harmonic"]){
                harmonicpotential = true;
                bonds[i].nominalLengthInNm=chromo.get_bond_nominal_length(i);
                bonds[i].stiffnessInKJPerNm2=chromo.get_spring_stiffness_coefficient();
            }
            
            
        }
        if(harmonicpotential){
            cout<<" Harmonic "<<endl;
            cout<<"\tCoeficient (KJ.Nm^-2.mol^-1 ) = " <<chromo.get_spring_stiffness_coefficient() <<endl;
        }
        vector<vector<int> > virtual_bond_list = chromo.get_virtual_bonds();
        
        int list_counter=0;
        for (int i=num_of_real_bonds+1; i<chromo_num_bonds; i++) {
            
            bonds[i].type = potentialModelIndex.Model["Virtual"];
            bonds[i].atoms[0]=virtual_bond_list[list_counter][0];
            bonds[i].atoms[1]=virtual_bond_list[list_counter][1];
            bonds[i].class_label = chromo.get_label() + chromo.get_label();
            bonds[i].nominalLengthInNm=0;
            bonds[i].stiffnessInKJPerNm2=0;
        }
    }
    return bonds;
}
