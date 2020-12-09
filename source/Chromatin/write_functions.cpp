//
//  write_functions.cpp
//  Mem
//
//  Created by Ali Farnudi on 19/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Chromatin.h"


using std::string;
using std::endl;
using std::cout;

void Chromatin::export_for_resume(int MD_step){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    write_resume_file<<num_of_node_types<<endl;
    for (int i=0; i<num_of_node_types; i++) {
        write_resume_file<<epsilon_LJ[i]<<"\t"<<sigma_LJ[i]<<endl;
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<ABC_index[i]<<"\n";
        //Node_force=0
    }
}

void Chromatin::export_for_resume(int MD_step, MyAtomInfo atoms[], int atom_count){
    std::ofstream write_resume_file;
    string resume_file_name="Results/Resumes/Resume_Chromatin_"+std::to_string(index)+"_";
    resume_file_name+=file_time;
    resume_file_name+=".txt";
    write_resume_file.open(resume_file_name.c_str());
    
    write_resume_file<<MD_step<<endl;
    write_resume_file<<Num_of_Nodes<<endl;
    
    write_resume_file<<num_of_node_types<<endl;
    for (int i=0; i<num_of_node_types; i++) {
        write_resume_file<<epsilon_LJ[i]<<"\t"<<sigma_LJ[i]<<endl;
    }
    
    for (int i=atom_count; i<atom_count+Num_of_Nodes; i++) {
        Node_Position[i-atom_count][0] = atoms[i].posInNm[0];
        Node_Position[i-atom_count][1] = atoms[i].posInNm[1];
        Node_Position[i-atom_count][2] = atoms[i].posInNm[2];
        
        Node_Velocity[i-atom_count][0] = atoms[i].velocityInNmperPs[0];
        Node_Velocity[i-atom_count][1] = atoms[i].velocityInNmperPs[1];
        Node_Velocity[i-atom_count][2] = atoms[i].velocityInNmperPs[2];
    }
    
    for (int i=0; i<Num_of_Nodes; i++) {
        write_resume_file<<Node_Position[i][0]<<"\t"<<Node_Position[i][1]<<"\t"<<Node_Position[i][2]<<"\n";
        write_resume_file<<Node_Velocity[i][0]<<"\t"<<Node_Velocity[i][1]<<"\t"<<Node_Velocity[i][2]<<"\n";
        write_resume_file<<ABC_index[i]<<"\n";
        //Node_force=0
    }
}

void Chromatin::export_coordinates(void){
    if (ExportGeneratedCoordinates) {
        std::string traj_name=generalParameters.trajectory_file_name+file_time+"_chromatin_"+std::to_string(index)+"_exported_coords.txt";
        
        std::ofstream exportcoords;
        exportcoords.open(traj_name.c_str());
        
        if (!exportcoords.is_open()){
            string errorMessage = TWARN;
            errorMessage+="I can't write the chromatin initial coordinates to storage. \nPath: "+traj_name;
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        for (int i=0; i<Num_of_Nodes; i++) {
            
            exportcoords<<std::fixed << std::setprecision(8)<<Node_Position[i][0]<<" "<<Node_Position[i][1]<<" "<<Node_Position[i][2]<<" "<<Node_Velocity[i][0]<<" "<<Node_Velocity[i][2]<<" "<<Node_Velocity[i][3]<<"\n";
        }
        exportcoords.close();
    }
    
}

