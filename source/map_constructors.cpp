#include "maps.hpp"
#include "General_constants.h"
using namespace std;

void read_general_parameters(map<string, double> general_param_map, string input_file_name){
    ifstream read_map(input_file_name.c_str());
    string param_name;
    double param_value;
    //    map<string, double> general_param_map;
    map<string, double>::iterator it;
    
    while (read_map>>param_name>>param_value) {
        general_param_map[param_name]=param_value;
    }
    
    param_name="MD_num_of_steps";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MD_num_of_steps=it->second;
    }
    param_name="MD_traj_save_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MD_traj_save_step=it->second;
    }
    param_name="MD_Time_Step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MD_Time_Step=it->second;
    }
    param_name="MD_KT";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MD_KT=it->second;
    }
    param_name="MD_thrmo_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MD_thrmo_step=it->second;
    }
    param_name="MC_step";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::MC_step=it->second;
    }
    param_name="Mem_fluidity";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::Mem_fluidity=it->second;
    }
    param_name="Lbox";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        GenConst::Lbox=it->second;
    }
    param_name="Periodic_condtion_status";
    it = general_param_map.find(param_name);
    if (it != general_param_map.end()){
        if (it->second==0.0) {
            GenConst::Periodic_condtion_status=false;
        } else {
            GenConst::Periodic_condtion_status=true;
        }
    }
}
 
