#include "OpenMM_funcs.hpp"

using std::set;
const int EndOfList=-1;

void set_bonded_forces(Bonds*                                 bonds,
                       OpenMM::HarmonicBondForce*            &HarmonicBond,
                       OpenMM::HarmonicBondForce*            &Kelvin_VoigtBond,
                       vector<OpenMM::CustomBondForce*>      &X4harmonics,
                       vector<OpenMM::CustomBondForce*>      &FENEs,
                       TimeDependantData*                    &time_dependant_data,
                       OpenMM::System                        &system
                       ){
    
    bool HarmonicBondForce=false;
    bool Kelvin_VoigtBondForce=false;
    
    set <std::string> FENE_classes;
    set <std::string> X4harmonic_classes;
    int FENE_index = -1;
    int X4harmonic_index = -1;
    
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
                
        }
        
        
    }
    
    if (HarmonicBondForce) {
        system.addForce(HarmonicBond);
    }
    
    if (Kelvin_VoigtBondForce) {
        system.addForce(Kelvin_VoigtBond);
        time_dependant_data->Kelvin_Voigt = true;
    }
}
