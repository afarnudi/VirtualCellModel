#include "maps.hpp"
#include "General_constants.h"

using namespace std;

void read_general_parameters(string input_file_name, vector<string> &membrane_config_list, vector<string> &chromatin_config_list){
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
        cout<<"Couldn't open the 'General_param_map.txt'. Please check the file.\n";
        exit(EXIT_FAILURE);
    }
    
    
    ifstream read_config_file(input_file_name.c_str());
    if (read_config_file.is_open()) {
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
//                        set_parameter(general_param_map, param_name, param_value);
                        //                general_param_map[param_name]=param_value;
                        for (int j=0; j<it->second; j++) {
                            cout<<split[i+2+j]<<endl<<endl;
                            membrane_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    } else if (it->first=="Num_of_Chromatins") {
                        //                        set_parameter(general_param_map, param_name, param_value);
                        //                general_param_map[param_name]=param_value;
                        for (int j=0; j<it->second; j++) {
                            cout<<split[i+2+j]<<endl<<endl;
                            chromatin_config_list.push_back(split[i+2+j]);
                        }
                        continue;
                    }
                    break;
                }
                //                cout<<split[i]<<"\t";
            }
            //            cout<<endl;
            
        }
//        for (auto const& x : general_param_map)
//        {
//            cout << x.first << " = " << x.second << endl;
//        }
//        while (read_config_file>>param_name>>param_value) {
//            if (param_name=="Num_of_Membranes") {
//                set_parameter(general_param_map, param_name, param_value);
////                general_param_map[param_name]=param_value;
//                for (int i=0; i<param_value; i++) {
//                    read_config_file>>param_name;
//                    membrane_list.push_back(param_name);
//                }
//                continue;
//            }
//            set_parameter(general_param_map, param_name, param_value);
////            general_param_map[param_name]=param_value;
//        }
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
    } else if (param_name=="MD_traj_save_step"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_traj_save_step=it->second;
        }
    } else if (param_name=="MD_Time_Step"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_Time_Step=it->second;
        }
    } else if (param_name=="MD_KT"){
        it = general_param_map.find(param_name);
        if (it != general_param_map.end()){
            GenConst::MD_KT=it->second;
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
    }
}
