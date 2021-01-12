#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

void set_bonded_forces(Bonds*                                 bonds,
                       OpenMM::HarmonicBondForce*            &HarmonicBond,
                       OpenMM::HarmonicBondForce*            &Kelvin_VoigtBond,
                       vector<OpenMM::CustomBondForce*>      &X4harmonics,
                       vector<OpenMM::CustomBondForce*>      &FENEs,
                       vector<OpenMM::CustomBondForce*>      &Gompperbond,
                       vector<OpenMM::CustomBondForce*>      &Gompperrep,
                       vector<OpenMM::CustomBondForce*>      &Contractiles,
                       vector<OpenMM::CustomBondForce*>      &KFs,
                       vector<OpenMM::CustomBondForce*>      &HillBonds,
                       vector<OpenMM::CustomBondForce*>      &Harmonic_minmax,
                       TimeDependantData*                    &time_dependant_data,
                       OpenMM::System                        &system
                       ){
    
    bool FENEBondForce=false;
    bool gompperBondForce=false;
    bool HarmonicBondForce=false;
    bool Kelvin_VoigtBondForce=false;
    bool Hill_Force=false;
    bool k_Force=false;
    
    set <std::string> FENE_classes;
    set <std::string> Gompper_classes;
    set <std::string> X4harmonic_classes;
    set <std::string> Contractile_classes;
    set <std::string> Hill_classes;
    set <std::string> KF_classes;
    set <std::string> Harmonic_minmax_classes;
    int FENE_index = -1;
    int Gompper_index = -1;
    int X4harmonic_index = -1;
    int Contractile_index = -1;
    int Hill_index = -1;
    int KF_index = -1;
    int Harmonic_minmax_index = -1;
    
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        
        const int*      atom = bonds[i].atoms;
        
        if (bonds[i].type == potentialModelIndex.Model["FENE"])
        {
            auto FENE_item = FENE_classes.find(bonds[i].class_label);
            if (FENE_item == FENE_classes.end()) {
                
                FENE_classes.insert(bonds[i].class_label);
                FENE_index++;
                
                FENEs.push_back(new OpenMM::CustomBondForce("-0.5*K_fene*R0*R0*log(1-(r*r/(R0*R0)))*step(0.99*R0-r)+4*epsilon_fene*((sigma_fene/r)^12-(sigma_fene/r)^6+0.25)*step(cut_fene-r)"));
                
                FENEs[FENE_index]->addPerBondParameter("K_fene");
                FENEs[FENE_index]->addPerBondParameter("R0");
                FENEs[FENE_index]->addPerBondParameter("epsilon_fene");
                FENEs[FENE_index]->addPerBondParameter("sigma_fene");
                FENEs[FENE_index]->addPerBondParameter("cut_fene");
                
                system.addForce(FENEs[FENE_index]);
            }
            vector<double> parameters={bonds[i].k_FENE_inKJpermol,
                                       bonds[i].FENER0inNm,
                                       bonds[i].epsilon_FENE_inKJpermol,
                                       bonds[i].nominalLengthInNm,
                                       bonds[i].nominalLengthInNm*(pow(2,1./6.))
                                       };
            
            FENEs[FENE_index]->addBond(atom[0], atom[1], parameters);
//            if (GenConst::Periodic_box) {
//                FENEs[FENE_index]->setUsesPeriodicBoundaryConditions(true);
//            }
        }
        else if (bonds[i].type == potentialModelIndex.Model["Gompper"])
        {
            auto Gompper_item = Gompper_classes.find(bonds[i].class_label);
            if (Gompper_item == Gompper_classes.end()) {
                
                Gompper_classes.insert(bonds[i].class_label);
                Gompper_index++;
                string label = bonds[i].class_label;

                string potentialbond = "step(0.99*"+label+"lmax-r)*"+"step(r-1.01*"+label+"lc0)*"+label+"bb*exp(1/("+label+"lc0-r))/("+label+"lmax-r)";

                Gompperbond.push_back(new OpenMM::CustomBondForce(potentialbond));
//                cout<<potentialbond<<endl;exit(0);
                Gompperbond[Gompper_index]->addGlobalParameter(label+"lc0" , bonds[i].gompperlc0);
                Gompperbond[Gompper_index]->addGlobalParameter(label+"bb"  , bonds[i].stiffnessInKJPerNm2);
                Gompperbond[Gompper_index]->addGlobalParameter(label+"lmax", bonds[i].gompperlmax);

                string potentialrep = "step(r-1.01*"+label+"lmin)*"+"step(0.99*"+label+"lc1-r)*"+label+"br*exp(1/(r-"+label+"lc1))/(r-"+label+"lmin)";
                Gompperrep.push_back(new OpenMM::CustomBondForce(potentialrep));
                cout<<potentialbond<<endl;cout<<potentialrep<<endl;
                Gompperrep[Gompper_index]->addGlobalParameter(label+"br"  , bonds[i].stiffnessInKJPerNm2);
                Gompperrep[Gompper_index]->addGlobalParameter(label+"lc1" , bonds[i].gompperlc1);
                Gompperrep[Gompper_index]->addGlobalParameter(label+"lmin", bonds[i].gompperlmin);
                    
//                string label = bonds[i].class_label;
//
//                string potentialbond = "step(r-1.01*lc0)*kgompper*exp(1/(lc0-r))/(lmax-r)";
//
//                Gompperbond.push_back(new OpenMM::CustomBondForce(potentialbond));
//
//                Gompperbond[Gompper_index]->addPerBondParameter("lc0");// , bonds[i].gompperlc0);
//                Gompperbond[Gompper_index]->addPerBondParameter("kgompper");//  , bonds[i].stiffnessInKJPerNm2);
//                Gompperbond[Gompper_index]->addPerBondParameter("lmax");//, bonds[i].gompperlmax);
//
//                string potentialrep = "step(0.99*lc1-r)*kgompper*exp(1/(r-lc1))/(r-lmin)";
//                Gompperrep.push_back(new OpenMM::CustomBondForce(potentialrep));
//                cout<<potentialbond<<endl;cout<<potentialrep<<endl;
//                Gompperrep[Gompper_index]->addPerBondParameter("kgompper");//  , bonds[i].stiffnessInKJPerNm2);
//                Gompperrep[Gompper_index]->addPerBondParameter("lc1");// , bonds[i].gompperlc1);
//                Gompperrep[Gompper_index]->addPerBondParameter("lmin");//, bonds[i].gompperlmin);
                
                
                system.addForce(Gompperbond[Gompper_index]);
//                system.addForce(Gompperrep [Gompper_index]);
            }
            
//            vector<double> bondparameters={bonds[i].gompperlc0,
//                                           bonds[i].stiffnessInKJPerNm2,
//                                           bonds[i].gompperlmax
//                                           };
//            Gompperbond[Gompper_index]->addBond(atom[0], atom[1], bondparameters);
            Gompperbond[Gompper_index]->addBond(atom[0], atom[1]);
            
//            vector<double> repparameters={bonds[i].stiffnessInKJPerNm2,
//                                          bonds[i].gompperlc1,
//                                          bonds[i].gompperlmin
//                                          };
//            Gompperrep [Gompper_index]->addBond(atom[0], atom[1], repparameters);
            Gompperrep [Gompper_index]->addBond(atom[0], atom[1]);

        }
        else if (bonds[i].type == potentialModelIndex.Model["Harmonic"])
        {
            HarmonicBondForce=true;
            // Note factor of 2 for stiffness below because Amber specifies the constant
            // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
            // it as used in the force term kx, with energy kx^2/2.
            HarmonicBond->addBond(atom[0], atom[1],
                                  bonds[i].nominalLengthInNm,
                                  bonds[i].stiffnessInKJPerNm2);
//            if (GenConst::Periodic_box) {
//                HarmonicBond->setUsesPeriodicBoundaryConditions(true);
//            }
        }
        else if (bonds[i].type == potentialModelIndex.Model["HarmonicX4"])
        {
            auto X4harmonic_item = X4harmonic_classes.find(bonds[i].class_label);
            if (X4harmonic_item == X4harmonic_classes.end()) {
                
                X4harmonic_classes.insert(bonds[i].class_label);
                X4harmonic_index++;
                
                X4harmonics.push_back( new OpenMM::CustomBondForce("k_bond*((r/r_rest)-1)^4"));
                
                X4harmonics[X4harmonic_index]->addPerBondParameter("r_rest");
                X4harmonics[X4harmonic_index]->addPerBondParameter("k_bond");
                system.addForce(X4harmonics[X4harmonic_index]);
            }
            double r_rest =bonds[i].nominalLengthInNm;
            double k_bond=bonds[i].stiffnessInKJPerNm4;
            vector<double> parameters;
            parameters.resize(2);
            parameters[0]=r_rest;
            parameters[1]=k_bond;
            X4harmonics[X4harmonic_index]->addBond(atom[0], atom[1], parameters);
//            if (GenConst::Periodic_box) {
//                X4harmonics[X4harmonic_index]->setUsesPeriodicBoundaryConditions(true);
//            }
        }
        else if (bonds[i].type == potentialModelIndex.Model["Kelvin-Voigt"])
        {
            Kelvin_VoigtBondForce=true;
            // Note factor of 2 for stiffness below because Amber specifies the constant
            // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
            // it as used in the force term kx, with energy kx^2/2.
            Kelvin_VoigtBond->addBond(atom[0], atom[1],
                                      bonds[i].nominalLengthInNm,
                                      bonds[i].stiffnessInKJPerNm2);
            
            time_dependant_data->Kelvin_Voigt_damp.push_back(bonds[i].dampInKJPsPerNm2);
//            if (GenConst::Periodic_box) {
//                Kelvin_VoigtBond->setUsesPeriodicBoundaryConditions(true);
//            }
        }
        else if (bonds[i].type == potentialModelIndex.Model["Contractile"])
        {
            auto Contractile_item = Contractile_classes.find(bonds[i].class_label);
            if (Contractile_item == Contractile_classes.end()) {
                
                Contractile_classes.insert(bonds[i].class_label);
                Contractile_index++;
                
                Contractiles.push_back( new OpenMM::CustomBondForce("(step(r-r_min) - step(r-r_max)) * F0 * (r - r0) + step(r-r_max) * F0 * (r_max - r0) + step(r_min - r) * F0 * (r_min - r0)"));
                
                Contractiles[Contractile_index]->addPerBondParameter("r0");
                Contractiles[Contractile_index]->addPerBondParameter("F0");
                Contractiles[Contractile_index]->addPerBondParameter("r_min");
                Contractiles[Contractile_index]->addPerBondParameter("r_max");
                system.addForce(Contractiles[Contractile_index]);
            }
            double r0 =bonds[i].nominalLengthInNm;
            double F0=bonds[i].F0;
            double r_min = bonds[i].r_min;
            double r_max = bonds[i].r_max;
            vector<double> parameters;
            parameters.resize(4);
            parameters[0]=r0;
            parameters[1]=F0;
            parameters[2]=r_min;
            parameters[3]=r_max;
            Contractiles[Contractile_index]->addBond(atom[0], atom[1], parameters);
            
//            if(GenConst::Periodic_box)
//            {
//                Contractiles[Contractile_index]->setUsesPeriodicBoundaryConditions(true);
//            }
            
            //time_dependant_data->hill_coefficient = bonds[i].hill_co;
            //time_dependant_data->Contractile_force = bonds[i].F0;
        }
        else if (bonds[i].type == potentialModelIndex.Model["Harmonic_minmax"])
        {
            auto Harmonic_minmax_item = Harmonic_minmax_classes.find(bonds[i].class_label);
            if (Harmonic_minmax_item == Harmonic_minmax_classes.end()) {
                
                Harmonic_minmax_classes.insert(bonds[i].class_label);
                Harmonic_minmax_index++;
                
                Harmonic_minmax.push_back( new OpenMM::CustomBondForce("(step(r-r_min) - step(r-r_max)) * (0.5 * k_bond *((r - r0)^2)) + step(r-r_max) * (0.5 * k_bond *((r_max - r0)^2)) + step(r_min - r) * (0.5 * k_bond *((r_min - r0)^2))"));
                
                Harmonic_minmax[Harmonic_minmax_index]->addPerBondParameter("r0");
                Harmonic_minmax[Harmonic_minmax_index]->addPerBondParameter("k_bond");
                Harmonic_minmax[Harmonic_minmax_index]->addPerBondParameter("r_min");
                Harmonic_minmax[Harmonic_minmax_index]->addPerBondParameter("r_max");
                system.addForce(Harmonic_minmax[Harmonic_minmax_index]);
            }
            double r0 =bonds[i].nominalLengthInNm;
            double k_bond=bonds[i].stiffnessInKJPerNm2;
            double r_min = bonds[i].r_min;
            double r_max = bonds[i].r_max;
            vector<double> parameters;
            parameters.resize(4);
            parameters[0]=r0;
            parameters[1]=k_bond;
            parameters[2]=r_min;
            parameters[3]=r_max;
            Harmonic_minmax[Harmonic_minmax_index]->addBond(atom[0], atom[1], parameters);
            
//            if(GenConst::Periodic_box)
//                      {
//                          Harmonic_minmax[Harmonic_minmax_index]->setUsesPeriodicBoundaryConditions(true);
//                      }
        } else if (bonds[i].type == potentialModelIndex.Model["hill"]){
            Hill_Force = true;
            auto Hill_item = Hill_classes.find(bonds[i].class_label);
            if (Hill_item == Hill_classes.end()) {
                
                Hill_classes.insert(bonds[i].class_label);
                Hill_index++;
                
                HillBonds.push_back( new OpenMM::CustomBondForce("(step(r-r_min) - step(r-r_max)) * F0 * (r - r0) + step(r-r_max) * F0 * (r_max - r0) + step(r_min - r) * F0 * (r_min - r0)"));
                
                HillBonds[Hill_index]->addPerBondParameter("r0");
                HillBonds[Hill_index]->addPerBondParameter("F0");
                HillBonds[Hill_index]->addPerBondParameter("r_min");
                HillBonds[Hill_index]->addPerBondParameter("r_max");
                HillBonds[Hill_index]->addPerBondParameter("b");
                system.addForce(HillBonds[Hill_index]);
            }
            double r0 =bonds[i].nominalLengthInNm;
            double F0=bonds[i].F0;
            double r_min = bonds[i].r_min;
            double r_max = bonds[i].r_max;
            double b = bonds[i].hill_co;
            vector<double> parameters;
            parameters.resize(5);
            parameters[0]=r0;
            parameters[1]=F0;
            parameters[2]=r_min;
            parameters[3]=r_max;
            parameters[4]=b;
            HillBonds[Hill_index]->addBond(atom[0], atom[1], parameters);
            
            time_dependant_data->hill_coefficient.push_back(b);
            time_dependant_data->hill_const_force.push_back(F0);
            
//            if(GenConst::Periodic_box)
//                      {
//                          HillBonds[Hill_index]->setUsesPeriodicBoundaryConditions(true);
//                      }
        }
        else if (bonds[i].type == potentialModelIndex.Model["KFs"]){
            k_Force = true;
            auto KF_item = KF_classes.find(bonds[i].class_label);
            if (KF_item == KF_classes.end()) {
                
                KF_classes.insert(bonds[i].class_label);
                KF_index++;
                
                KFs.push_back( new OpenMM::CustomBondForce("(step(r-r_min) - step(r-r_max)) * F0 * (r - r0) + step(r-r_max) * F0 * (r_max - r0) + step(r_min - r) * F0 * (r_min - r0)"));
                
                KFs[KF_index]->addPerBondParameter("r0");
                KFs[KF_index]->addPerBondParameter("F0");
                KFs[KF_index]->addPerBondParameter("r_min");
                KFs[KF_index]->addPerBondParameter("r_max");
                KFs[KF_index]->addPerBondParameter("n");
                KFs[KF_index]->addPerBondParameter("Fc");
                system.addForce(KFs[KF_index]);
            }
            double r0 =bonds[i].nominalLengthInNm;
            double F0=bonds[i].F0;
            double r_min = bonds[i].r_min;
            double r_max = bonds[i].r_max;
            double n = bonds[i].k_F0;
            double Fc=bonds[i].F0;
            vector<double> parameters;
            parameters.resize(6);
            parameters[0]=r0;
            parameters[1]=F0;
            parameters[2]=r_min;
            parameters[3]=r_max;
            parameters[4]=n;
            parameters[5]=Fc;
            KFs[KF_index]->addBond(atom[0], atom[1], parameters);
            
            //time_dependant_data->hill_coefficient.push_back(b);
            //time_dependant_data->hill_const_force.push_back(F0);
        }
        
        
    }
    
    if (HarmonicBondForce) {
        system.addForce(HarmonicBond);
    }
    
    if (Kelvin_VoigtBondForce) {
        system.addForce(Kelvin_VoigtBond);
        time_dependant_data->Kelvin_Voigt = true;
    }
    
    if (Hill_Force) {
        time_dependant_data->HillForce = true;
    }
    
    if (k_Force) {
        time_dependant_data->kForce = true;
    }
}
