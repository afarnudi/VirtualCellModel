#include "maps.hpp"
#include "General_constants.h"

using namespace std;

void read_general_parameters(string input_file_name, vector<string> &membrane_config_list, vector<string> &chromatin_config_list, vector<string> &actin_config_list, vector<string> &ecm_config_list, vector<string> &pointparticle_config_list){
    ifstream read_map("General_param_map.txt");
    map<string, double> general_param_map;
    map<string, double>::iterator it;
    string param_name;
    double param_value;
    
    if (read_map.is_open()) {
        //    map<string, double> general_param_map;
        
        while (read_map>>param_name>>param_value) {
            general_param_map[param_name]=param_value;
            set_parameter(general_param_map, param_name, param_value);
        }
    } else {
        cout<<"Couldn't open the 'General_param_map.txt'. Please check the file and make sure that the file is in the same directory as the executable file.\n";
        exit(EXIT_FAILURE);
    }
    
    
    ifstream read_config_file(input_file_name.c_str());
    if (read_config_file.is_open()) {
        if (GenConst::Testmode) {
            cout<<"\nGeneral Parameter file opened successfully.\nList of configuration files:\n";
        }
        
        string line;
        int line_num=0;
        string comment="//";
        //        char delimiter=' ';
        while(getline(read_config_file, line)){
            line_num++;
            
            if(line.empty()){
                continue;
            }
            
            istringstream iss(line);
            vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
            
            for (int i=0; i<split.size(); i++) {
                if (split[i] == comment || (split[i][0]=='/' && split[i][1]=='/')) {
                    break;
                }
                
                it = general_param_map.find(split[i]);
                if (it != general_param_map.end()){
                    it->second=stod(split[i+1]);
                    set_parameter(general_param_map, it->first, it->second);
                    if (it->first=="Num_of_Membranes") {
                        
                        for (int j=0; j<it->second; j++) {
                            if (GenConst::Testmode) {
                                cout<<"\t"<<split[i+2+j]<<endl;
                            }
                            
                            membrane_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="Num_of_Chromatins") {
                        
                        for (int j=0; j<it->second; j++) {
                            cout<<"\t"<<split[i+2+j]<<endl;
                            chromatin_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="Num_of_Actins") {
                        
                        for (int j=0; j<it->second; j++) {
                            cout<<"\t"<<split[i+2+j]<<endl;
                            actin_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="Num_of_ECMs") {
                        
                        for (int j=0; j<it->second; j++) {
                            cout<<"\t"<<split[i+2+j]<<endl;
                            ecm_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="Num_of_pointparticles") {
                        
                        for (int j=0; j<it->second; j++) {
                            cout<<"\t"<<split[i+2+j]<<endl;
                            pointparticle_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="trajectory_file_name") {
                        if (it->second!=0) {
                            GenConst::trajectory_file_name=split[i+2];
                        } else {
                            GenConst::trajectory_file_name="VCProject_";
                        }
                        
                        continue;
                    } else if (it->first=="Interaction_map") {
                        if (it->second==0) {
                            GenConst::Interaction_map = false;
                            GenConst::Interaction_map_file_name="interaction_map.txt";
                        } else {
                            GenConst::Interaction_map = true;
                            GenConst::Interaction_map_file_name=split[i+2];
                        }
                        
                        continue;
                    } else if (it->first=="Membrane_label") {
                        if (it->second==0) {
                            GenConst::Membrane_label="mem";
                        } else {
                            GenConst::Membrane_label=split[i+2];
                        }
                        
                        continue;
                    } else if (it->first=="Actin_label") {
                        if (it->second==0) {
                            GenConst::Actin_label="act";
                        } else {
                            GenConst::Actin_label=split[i+2];
                        }
                        
                        continue;
                    } else if (it->first=="Chromatin_label") {
                        if (it->second==0) {
                            GenConst::Chromatin_label="chr";
                        } else {
                            GenConst::Chromatin_label=split[i+2];
                        }
                        
                        continue;
                    } else if (it->first=="ECM_label") {
                        if (it->second==0) {
                            GenConst::ECM_label="ecm";
                        } else {
                            GenConst::ECM_label=split[i+2];
                        }
                        
                        continue;
                    } else if (it->first=="Load_from_checkpoint") {
                        if (it->second==0) {
                            GenConst::Load_from_checkpoint = false;
                            GenConst::Checkpoint_file_name = "None";
                        } else {
                            GenConst::Load_from_checkpoint = true;
                            GenConst::Checkpoint_file_name=split[i+2];
                        }
                        
                        continue;
                    }  else if (it->first=="Checkpoint_path") {
                        if (it->second==0) {
                            GenConst::Checkpoint_path = "Results/Resumes/OpenMM/";
                        } else {
                            GenConst::Checkpoint_path = split[i+2];
                        }
                        
                        continue;
                    }
                    
                    break;
                }
                //                cout<<split[i]<<"\t";
            } // End of for (int i=0; i<split.size(); i++) {
            //            cout<<endl;
        } //End of while(getline(read_config_file, line)){
        if (GenConst::Report_Interval_In_Fs==0) {
            GenConst::Report_Interval_In_Fs = GenConst::MD_traj_save_step*GenConst::Step_Size_In_Fs;
        }
        if (GenConst::Simulation_Time_In_Ps==0) {
            GenConst::Simulation_Time_In_Ps = GenConst::MD_num_of_steps*GenConst::Step_Size_In_Fs;
        }
//        cout<<endl;
    } else {
        cout<<"Couldn't open the config file.\n";
        exit(EXIT_FAILURE);
    }
    
}

void set_parameter(map<string, double> &general_param_map, string param_name, double param_value){
    
    map<string, double>::iterator it;
    if (param_name=="MD_num_of_steps") {
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_num_of_steps=it->second;
        }
    } else if (param_name=="Simulation_Time_In_Ps"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Simulation_Time_In_Ps=it->second;
        }
    } else if (param_name=="MD_traj_save_step"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_traj_save_step=it->second;
        }
    } else if (param_name=="Report_Interval_In_Fs"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Report_Interval_In_Fs=it->second;
        }
    } else if (param_name=="Step_Size_In_Fs"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Step_Size_In_Fs=it->second;
        }
    } else if (param_name=="MD_T"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_T=it->second;
            GenConst::Buffer_temperature=GenConst::MD_T;
        }
    } else if (param_name=="K"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::K=it->second;
        }
    } else if(param_name=="MD_thrmo_step"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_thrmo_step=it->second;
        }
    } else if (param_name=="MC_step"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MC_step=it->second;
        }
    } else if (param_name=="Mem_fluidity"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Mem_fluidity=it->second;
        }
    } else if(param_name=="Lbox"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Lbox=it->second;
            if (it->second<=0.1) {
                GenConst::Periodic_box=false;
            } else {
                GenConst::Periodic_box=true;
            }
        }
    } else if (param_name=="Num_of_Membranes"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Num_of_Membranes=it->second;
        } else {
            GenConst::Num_of_Membranes=0;
        }
    } else if (param_name=="Num_of_Chromatins"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Num_of_Chromatins=it->second;
        } else {
            GenConst::Num_of_Chromatins=0;
        }
    } else if (param_name=="Num_of_Actins"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Num_of_Actins=it->second;
        } else {
            GenConst::Num_of_Actins=0;
        }
    } else if (param_name=="Num_of_ECMs"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Num_of_ECMs=it->second;
        } else {
            GenConst::Num_of_ECMs=0;
        }
    }else if (param_name=="Num_of_pointparticles"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Num_of_pointparticles=it->second;
        } else {
            GenConst::Num_of_pointparticles=0;
        }
    } else if (param_name=="Bussi_tau"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Bussi_tau=it->second;
        }
    } else if (param_name=="Actin_Membrane_Bond_Coefficient"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Actin_Membrane_Bond_Coefficient=it->second;
        }
    } else if (param_name=="excluded_volume_interaction"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::Excluded_volume_interaction= false;
            } else {
                GenConst::Excluded_volume_interaction= true;
            }
        }
    } else if (param_name=="Interaction_map"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::Interaction_map= false;
            } else {
                GenConst::Interaction_map= true;
            }
        }
    } else if (param_name=="Membrane_label"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Membrane_label="mem";
        }
    } else if (param_name=="Actin_label"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Actin_label="act";
        }
    } else if (param_name=="Chromatin_label"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Chromatin_label="chr";
        }
    } else if (param_name=="ECM_label"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::ECM_label="ecm";
        }
    } else if (param_name=="OpenMM"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::OpenMM= false;
            } else {
                GenConst::OpenMM= true;
            }
        }
    } else if (param_name=="sigma_LJ_12_6"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::sigma_LJ_12_6=it->second;
        }
    } else if (param_name=="epsilon_LJ_12_6"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::epsilon_LJ_12_6=it->second;
        }
    } else if (param_name=="Integrator_type"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Integrator_type=it->second;
        }
    } else if (param_name=="frictionInPs"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::frictionInPs=it->second;
        }
    } else if (param_name=="temperature"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::temperature=it->second;
        }
    } else if (param_name=="WantEnergy"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::WantEnergy= false;
            } else {
                GenConst::WantEnergy= true;
            }
            
        }
    } else if (param_name=="WantForce"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::WantForce= false;
            } else {
                GenConst::WantForce= true;
            }
            
        }
    } else if (param_name=="WriteVelocitiesandForces"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::WriteVelocitiesandForces= false;
            } else {
                GenConst::WriteVelocitiesandForces= true;
            }
        }
    } else if (param_name=="Load_from_checkpoint"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Load_from_checkpoint=it->second;
        }
    } else if (param_name=="Checkpoint_path"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::Checkpoint_path="Results/Resumes/OpenMM/";
        }
    } else if (param_name=="write_bonds_to_PDB"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::write_bonds_to_PDB= false;
            } else {
                GenConst::write_bonds_to_PDB= true;
            }
        }
    } else if (param_name=="CMMotionRemover"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::CMMotionRemover= false;
            } else {
                GenConst::CMMotionRemover= true;
            }
        }
    } else if (param_name=="CMMotionRemoverStep"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::CMMotionRemoverStep=it->second;
        }
    } else if (param_name=="CreateCheckpoint"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::CreateCheckpoint = false;
            } else {
                GenConst::CreateCheckpoint = true;
            }
        }
    } else if (param_name=="Wantvoronoi"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::Wantvoronoi = false;
            } else {
                GenConst::Wantvoronoi = true;
            }
        }
    } else if (param_name=="Testmode"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second == 0) {
                GenConst::Testmode = false;
            } else {
                GenConst::Testmode = true;
            }
        }
    } else if (param_name=="MCBarostatPressure"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MCBarostatPressure = it->second;
        }
    } else if (param_name=="MCBarostatTemperature"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MCBarostatTemperature = it->second;
        }
    } else if (param_name=="MCBarostatFrequency"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MCBarostatFrequency = it->second;
        }
    }
    
    
}
