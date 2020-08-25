#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

void set_bonded_forces(Bonds*                                 bonds,
                       OpenMM::HarmonicBondForce*            &HarmonicBond,
                       OpenMM::HarmonicBondForce*            &Kelvin_VoigtBond,
                       vector<OpenMM::CustomBondForce*>      &X4harmonics,
                       vector<OpenMM::CustomBondForce*>      &FENEs,
                       vector<OpenMM::CustomBondForce*>      &Contractiles,
                       vector<OpenMM::CustomBondForce*>      &KFs,
                       vector<OpenMM::CustomBondForce*>      &HillBonds,
                       vector<OpenMM::CustomBondForce*>      &Harmonic_minmax,
                       TimeDependantData*                    &time_dependant_data,
                       OpenMM::System                        &system
                       ){
    
    bool HarmonicBondForce=false;
    bool Kelvin_VoigtBondForce=false;
    bool Hill_Force=false;
    bool k_Force=false;
    
    set <std::string> FENE_classes;
    set <std::string> X4harmonic_classes;
    set <std::string> Contractile_classes;
    set <std::string> Hill_classes;
    set <std::string> KF_classes;
    set <std::string> Harmonic_minmax_classes;
    int FENE_index = -1;
    int X4harmonic_index = -1;
    int Contractile_index = -1;
    int Hill_index = -1;
    int KF_index = -1;
    int Harmonic_minmax_index = -1;
    
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        
        const int*      atom = bonds[i].atoms;
        
        switch (bonds[i].type) {
            case 1://FENE
            {
                
                
                auto FENE_item = FENE_classes.find(bonds[i].class_label);
                if (FENE_item == FENE_classes.end()) {
                    
                    FENE_classes.insert(bonds[i].class_label);
                    FENE_index++;
                    
                    FENEs.push_back(new OpenMM::CustomBondForce("4*ep*( (s/r)^12-(s/r)^6 + 0.25) -0.5*K_fene*R*R*log(1-(r*r/(R*R) ))"));
                    
                    FENEs[FENE_index]->addPerBondParameter("s");
                    FENEs[FENE_index]->addPerBondParameter("R");
                    FENEs[FENE_index]->addPerBondParameter("ep");
                    FENEs[FENE_index]->addPerBondParameter("K_fene");
                    
                    system.addForce(FENEs[FENE_index]);
                }
                vector<double> parameters={bonds[i].FENE_lmininNm,
                                           bonds[i].FENE_lmaxinNm,
                                           bonds[i].FENE_epsilon,
                                           bonds[i].FENE_k
                                           };
                
                FENEs[FENE_index]->addBond(atom[0], atom[1], parameters);
                if (GenConst::Periodic_box) {
                    FENEs[FENE_index]->setUsesPeriodicBoundaryConditions(true);
                }
            }
                break;
            case 2://Harmonic
            {
                HarmonicBondForce=true;
                // Note factor of 2 for stiffness below because Amber specifies the constant
                // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
                // it as used in the force term kx, with energy kx^2/2.
                HarmonicBond->addBond(atom[0], atom[1],
                                      bonds[i].nominalLengthInNm,
                                      bonds[i].stiffnessInKJPerNm2);
                if (GenConst::Periodic_box) {
                    HarmonicBond->setUsesPeriodicBoundaryConditions(true);
                }
                
                
            }
                break;
            case 3:// X4harmonic
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
                if (GenConst::Periodic_box) {
                    X4harmonics[X4harmonic_index]->setUsesPeriodicBoundaryConditions(true);
                }
            }
                
                
                break;
                
                
            case 4://Kelvin-Voigt
            {
                Kelvin_VoigtBondForce=true;
                // Note factor of 2 for stiffness below because Amber specifies the constant
                // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
                // it as used in the force term kx, with energy kx^2/2.
                Kelvin_VoigtBond->addBond(atom[0], atom[1],
                                          bonds[i].nominalLengthInNm,
                                          bonds[i].stiffnessInKJPerNm2);
                
                time_dependant_data->Kelvin_Voigt_damp.push_back(bonds[i].dampInKJPsPerNm2);
                if (GenConst::Periodic_box) {
                    Kelvin_VoigtBond->setUsesPeriodicBoundaryConditions(true);
                }
            }
                break;
                
                
            case 6:// Contractile constant
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
                
                if(GenConst::Periodic_box)
                {
                    Contractiles[Contractile_index]->setUsesPeriodicBoundaryConditions(true);
                }
                
                //time_dependant_data->hill_coefficient = bonds[i].hill_co;
                //time_dependant_data->Contractile_force = bonds[i].F0;
            }
                
                
                break;
                
                
            case 7:// Harmonic_minmax
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
                
                if(GenConst::Periodic_box)
                          {
                              Harmonic_minmax[Harmonic_minmax_index]->setUsesPeriodicBoundaryConditions(true);
                          }
            }
                
                
                break;
                
                
                
                case 8:// hill force
                {
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
                    
                    if(GenConst::Periodic_box)
                              {
                                  HillBonds[Hill_index]->setUsesPeriodicBoundaryConditions(true);
                              }
                }
                    
                    
                    break;
                
                
                
                case 9:// KFs
                {
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
                    
                    
                    break;

                
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
