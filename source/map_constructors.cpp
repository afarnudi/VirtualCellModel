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
        cout<<"General Parameter file opened successfully.\nList of config-file:\n";
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
                            cout<<"\t"<<split[i+2+j]<<endl;
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
//                            cout<<"\t"<<split[i+2+j]<<endl;
                            ecm_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    }  else if (it->first=="Num_of_pointparticles") {

                        for (int j=0; j<it->second; j++) {
                            cout<<"\t"<<split[i+2+j]<<endl;
                            pointparticle_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    }else if (it->first=="trajectory_file_name") {
                        if (it->second!=0) {
                            GenConst::trajectory_file_name=split[i+2];
                        } else {
                            GenConst::trajectory_file_name="VCProject_";
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
        cout<<endl;
    } else {
        cout<<"Couldn't open the config file.\n";
        exit(EXIT_FAILURE);
    }

//    cout<<"Mem_fluidity "<<GenConst::Mem_fluidity<<"\nMD_num_of_steps "<<GenConst::MD_num_of_steps<<"\nNum_of_Membranes "<<GenConst::Num_of_Membranes<<"\nMD_traj_save_step "<<GenConst::MD_traj_save_step<<"\nMembrane file name "<<membrane_list[0]<<endl;
//    exit(EXIT_SUCCESS);

//    for (auto const& x : general_param_map)
//    {
//        cout << x.first << " = " << x.second << endl;
//    }
//    for (int i=0; i<membrane_list.size(); i++) {
//        cout<<membrane_list[i]<<"\n";
//    }
//    exit(EXIT_SUCCESS);
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
//            cout<<"Simulation_Time_In_Ps = "<< GenConst::Simulation_Time_In_Ps<<endl;
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
        }
    } else if(param_name=="Periodic_condtion_status"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            if (it->second==0.0) {
                GenConst::Periodic_condtion_status=false;
            } else {
                GenConst::Periodic_condtion_status=true;
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
    }

}
