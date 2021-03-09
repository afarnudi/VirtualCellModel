#include "Membrane.h"
#include "General_functions.hpp"
#include "Global_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"
#include <vector>
#include <stdlib.h>
#include <math.h>


using OpenMM::Vec3;
//                           PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
void myWritePDBFrame(int frameNum,
                     double timeInPs,
                     double energyInKJ,
                     double potential_energyInKJ,
                     const MyAtomInfo atoms[],
                     const Bonds bonds[])
{
    int EndOfList=-1;
    
    string traj_name= generalParameters.trajectory_file_name+".pdb";
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", frameNum);
    fprintf(pFile,"REMARK 250 time=%.3f ps; energy=%6.6f potential energy=%.3f KJ/mole\n",
            timeInPs,
            energyInKJ,
            potential_energyInKJ);
    //    cout<<endl;
    //    cout<<"kbend*angle "<<potential_energyInKJ<<endl;
    //    cout<<"thetha "<<potential_energyInKJ*180/M_PI<<endl;
    int index=0;
    string hist = atoms[0].pdb;
    if (atoms[0].class_label == "Chromatin") {
        hist.pop_back();
    }
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    //    double occ=1;
    for (int n=0; atoms[n].type != EndOfList; ++n){
        string new_label = atoms[n].pdb;
        if (atoms[n].class_label == "Chromatin") {
            new_label.pop_back();
        }
        if (atoms[n].class_label == "ECM") {
            new_label.pop_back();
        }
        
        if (hist != new_label) {
            index++;
            hist = new_label;
        }
        //        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
        fprintf(pFile,"ATOM %6d %4s ETH %c%4.0f    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                atoms[n].pdb,
                chain[index],
                double(index),
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2],
                atoms[n].stretching_energy,
                atoms[n].energyInKJ);//,
        //                atoms[n].symbol);
    }
    
    
    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}

void writeXYZFrame  (int atom_count,
                     const MyAtomInfo atoms[],
                     std::string traj_name)
{
    int EndOfList=-1;
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"%d\ncomment\n", atom_count);
    
    for (int n=0; atoms[n].type != EndOfList; ++n){
        fprintf(pFile,"%4s\t%8.3f\t%8.3f\t%8.3f\n",
                atoms[n].pdb,
                atoms[n].posInNm[0],
                atoms[n].posInNm[1],
                atoms[n].posInNm[2]);
    }
    fclose (pFile);
}







void myWritePSF(int   num_of_atoms,
                int   num_of_bonds,
                const MyAtomInfo   atoms[],
                const Bonds        bonds[])
{
    int EndOfList=-1;
    
    string traj_name= generalParameters.trajectory_file_name+".psf";
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile," PSF\n\n");
    fprintf(pFile,"       1  !TITLE\n");
    fprintf(pFile," vmd files (psf,pdb) for VCM\n\n");
    fprintf(pFile,"   %5d !NATOM\n",num_of_atoms);
    
    
    int index=0;
    string hist = atoms[0].pdb;
    if (atoms[0].class_label == "Chromatin") {
        hist.pop_back();
    }
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    
    //    double occ=1;
    for (int n=0; atoms[n].type != EndOfList; ++n){
        string new_label = atoms[n].pdb;
        if (atoms[n].class_label == "Chromatin") {
            new_label.pop_back();
        }
        if (atoms[n].class_label == "ECM") {
            new_label.pop_back();
        }
        
        if (hist != new_label) {
            index++;
            hist = new_label;
        }
        //        fprintf(pFile,"ATOM  %5d %4s ETH %c   %4.0f %8.3f%8.3f%8.3f%6.2f%6.1f          %c\n",
        fprintf(pFile,"   %5d POLY%5d POLY %4s CEL1 %10.2f    %10.2f 1\n",
                n+1,
                index,
                atoms[n].pdb,
                1.0,
                atoms[n].mass);
    }
    num_of_bonds=0;
    for (int n=0; bonds[n].type != EndOfList; ++n){
        if (bonds[n].type != potentialModelIndex.Model["None"]) {
            num_of_bonds++;
        }
        
    }
    
    
    fprintf(pFile,"\n");
    fprintf(pFile,"   %5d !NBOND\n",num_of_bonds);
    int endline_counter=0;
    for (int n=0; bonds[n].type != EndOfList; ++n){
        if (bonds[n].type != potentialModelIndex.Model["None"]) {
            fprintf(pFile,"   %5d   %5d",
                    bonds[n].atoms[0]+1,
                    bonds[n].atoms[1]+1);
            endline_counter++;
            if (endline_counter>3) {
                endline_counter=0;
                fprintf(pFile,"\n");
            }
        }
        
    }
    fclose (pFile);
}
