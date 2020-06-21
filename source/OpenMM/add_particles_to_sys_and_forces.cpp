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
        
        const Vec3 posInNm(atoms[n].initPosInNm[0],
                           atoms[n].initPosInNm[1],
                           atoms[n].initPosInNm[2]);
        const Vec3 velocityInNmperPs(atoms[n].velocityInNmperPs[0],
                                     atoms[n].velocityInNmperPs[1],
                                     atoms[n].velocityInNmperPs[2]);
        
        initialPosInNm.push_back(posInNm);
        initialVelInNmperPs.push_back(velocityInNmperPs);
        
        //add particles to the excluded volume force. The number of particles should be equal to the number particles in the system. The exluded interaction lists should be defined afterwards.
        std::vector<double> sigma_ev;
        sigma_ev.push_back(atoms[n].radius);
        
        for (int i=0; i<ExcludedVolumes.size(); i++) {
            ExcludedVolumes[i]->addParticle();
        }
        for (int i=0; i<LJ_12_6_interactions.size(); i++) {
            LJ_12_6_interactions[i]->addParticle();
        }
        
        
        
    }
}
