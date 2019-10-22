#include "OpenMM_funcs.hpp"


const int EndOfList=-1;
using OpenMM::Vec3;
void add_particles_to_system_and_forces(const MyAtomInfo                       atoms[],
                                        vector<Vec3>                          &initialPosInNm,
                                        vector<Vec3>                          &initialVelInNmperPs,
                                        vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                                        vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                                        OpenMM::System                        &system){
    int EndOfList=-1;
    std::vector<double> sigma_ev;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        //        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atoms[n].mass);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        const Vec3 velocityInAngperPs(atoms[n].velocityInAngperPs[0],
                                      atoms[n].velocityInAngperPs[1],
                                      atoms[n].velocityInAngperPs[2]);
        //        cout<<atoms[n].velocityInAngperPs[0]<<"\t"<<
        //              atoms[n].velocityInAngperPs[1]<<"\t"<<
        //              atoms[n].velocityInAngperPs[2]<<endl;
        initialPosInNm.push_back(posInNm);
        initialVelInNmperPs.push_back(velocityInAngperPs);
        
        //add particles to the excluded volume force. The number of particles should be equal to the number particles in the system. The exluded interaction lists should be defined afterwards.
        std::vector<double> sigma_ev;
        sigma_ev.push_back(atoms[n].radius
                           * OpenMM::NmPerAngstrom);
        for (int i=0; i<ExcludedVolumes.size(); i++) {
            ExcludedVolumes[i]->addParticle(sigma_ev);
        }
        for (int i=0; i<LJ_12_6_interactions.size(); i++) {
            LJ_12_6_interactions[i]->addParticle();
        }
        
        
        
    }
}
