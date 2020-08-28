#include <sstream>
#include "Chromatin.h"

using namespace std;

void Chromatin::import_config(string config_file_name, double min_radius){
    
    map<string, double>::iterator it;
    string resume_file_name;
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    
    if (read_config_file.is_open()) {
        cout<<"'"<<TFILE<<config_file_name<<TRESET<<"' file opened "<<TSUCCESS<<"successfully"<<TRESET<<".";
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
            
            if (split[0] == comment || (split[0][0]=='/' && split[0][1]=='/')) {
                continue;
            }
            
            param_map[split[0]]=stod(split[1]);
            
            if (split[0]=="Resume") {
                
                if (stoi(split[1])==0) {
                    cout<<"Resume flag off. The Chromatins will be initiated using the config parameters.\n";
                } else {
                    resume=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag on. Chromatin will resume using the '"<<resume_file_name<<" file.\n";
                }
            } else {
                set_map_parameter(split[0], param_map[split[0]]);
                
            }
            
        }//End of while(getline(read_config_file, line))
        
        
    } else {//End of if (read_config_file.is_open())
        cout<<"Couldn't open the '"<<config_file_name<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume) {
        import_resume(resume_file_name);
    } else {
        if (Num_of_Nodes==0) {
            cout<< "Error. \nPlease specify the number of Chromatin nodes in the "<<config_file_name<< " file.\n";
            exit(EXIT_FAILURE);
        }
        initialise(min_radius);
    }
    //        else {
    //            cout<<"Resume is off and no meshfile name is provided for initilisation. Please check the membrane config file.\n";
    //        }
}

void Chromatin::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    string resume_file_name;
    string import_file_name;
    ifstream read_config_file(config_file_name.c_str());
    bool resume_flag=false;
    bool import_flag=false;
    GenConst::ChromatinVirtualSites = false;
    
    if (read_config_file.is_open()) {
        cout<<"'"<<config_file_name<<"' file opened successfully.\n";
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
            
            if (split[0] == comment || (split[0][0]=='/' && split[0][1]=='/')) {
                continue;
            }
            
            param_map[split[0]]=stod(split[1]);
            
            if (split[0]=="Resume") {
                
                if (stoi(split[1])==0) {
                    cout<<"Resume flag "<<TRED<<"off"<<TRESET<<".\n";
                } else {
                    resume_flag=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag on. Chromatin will resume using the '"<<resume_file_name<<" file.\n";
                }
            } else if(split[0]=="epsilon"){
                if (split.size()<num_of_node_types + 1) {
                    cout<<"Too few arguments for the Lenard Jones epsilon interaction of node types.\nNeed "<<num_of_node_types<<" argumens.\n "<<split.size()-1<<" was provided.\nWill resume chromatin interactions with default value, 0.\n";
                    
                    for (int i=0; i<split.size(); i++) {
                        epsilon_LJ[i]=stod(split[i+1]);
                    }
                } else {
                    for (int i=0; i<num_of_node_types; i++) {
                        epsilon_LJ[i]=stod(split[i+1]);
                    }
                }
            } else if(split[0]=="sigma"){
                if (split.size()<num_of_node_types + 1) {
                    cout<<"Too few arguments for the Lenard Jones sigma interaction of node types.\nNeed "<<num_of_node_types<<" argumens.\n "<<split.size()-1<<" was provided.\nWill resume chromatin interactions with default value, 1.5x node_radius \n";
                    for (int i=0; i<split.size(); i++) {
                        sigma_LJ[i]=stod(split[i+1]);
                    }
                } else {
                    for (int i=0; i<num_of_node_types; i++) {
                        sigma_LJ[i]=stod(split[i+1]);
                    }
                }
            } else if(split[0]=="Import_coordinates"){
                if (stoi(split[1])==0) {
                    cout<<"Import flag off. The Chromatins will be initiated using the config parameters.\n";
                } else {
                    import_flag=true;
                    import_file_name=split[2];
                    cout<<"Import flag "<<TGREEN<<"on"<<TRESET<<". importing coordinatesfrom '"<<TFILE<<import_file_name<<TRESET<<".\n";
                }
            } else {
                set_map_parameter(split[0], param_map[split[0]]);
                
            }
            
        }//End of while(getline(read_config_file, line))
        
        if (epsilon_LJ.size() != num_of_node_types and num_of_node_types>1) {
            cout<<"error: You need so specify an epsilon interaction value for chromatin chains with more than one node type. (set in the chromatin configuration file)\n";
            exit(EXIT_FAILURE);
        }
        
        if (sigma_LJ.size() != num_of_node_types and num_of_node_types>1) {
            cout<<"All Lenard Jones 12 6 sigmas set to default value (1.5 x Node Radius)\n";
            sigma_LJ.resize(num_of_node_types, 1.5*Node_radius);
        }
        
        if (num_of_node_types == 1) {
            if (epsilon_LJ.size() == 0) {
                epsilon_LJ.resize(1,0);
            }
            if (sigma_LJ.size() == 0) {
                sigma_LJ.resize(1, 1.5*Node_radius);
            }
        }
        
    } else {//End of if (read_config_file.is_open())
        cout<<"Couldn't open the '"<<config_file_name<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume_flag) {
        import_resume(resume_file_name);
    } else if (import_flag) {
        import_coordinates(import_file_name);
    } else {
        if (Num_of_Nodes==0) {
            cout<< "Error. \nPlease specify the number of Chromatin nodes in the "<<config_file_name<< " file.\n";
            exit(EXIT_FAILURE);
        }
        initialise();
    }
    //        else {
    //            cout<<"Resume is off and no meshfile name is provided for initilisation. Please check the membrane config file.\n";
    //        }
}


void Chromatin::set_map_parameter(string param_name, double param_value){
    
    //    map<string, double>::iterator it;
    if (param_name=="Node_Mass") {
        Node_Mass=param_value;
    } else if (param_name=="Node_radius"){
        Node_radius=param_value;
    } else if (param_name=="spring_model"){
        spring_model=param_value;
    } else if (param_name=="Spring_coefficient"){
        Spring_coefficient=param_value;
    } else if (param_name=="Damping_coefficient"){
        Damping_coefficient=param_value;
    } else if (param_name=="Shift_in_X_direction"){
        Shift_in_X_direction=param_value;
    } else if (param_name=="Shift_in_Y_direction"){
        Shift_in_Y_direction=param_value;
    } else if (param_name=="Shift_in_Z_direction"){
        Shift_in_Z_direction=param_value;
    } else if (param_name=="x_speed"){
        x_speed=param_value;
    } else if (param_name=="y_speed"){
        y_speed=param_value;
    } else if (param_name=="z_speed"){
        z_speed=param_value;
    } else if (param_name=="rescale_factor"){
        rescale_factor=param_value;
    } else if (param_name=="Num_of_Nodes" && Num_of_Nodes==0){
        Num_of_Nodes=param_value;
    } else if (param_name=="num_of_node_types"){
        num_of_node_types=param_value;
        epsilon_LJ.resize(num_of_node_types,0);
        sigma_LJ.resize(num_of_node_types,2.5*Node_radius);
    } else if (param_name=="bond_length"){
        bond_length=param_value;
    } else if (param_name=="bond_radius"){
        bond_radius=param_value;
    } else if (param_name=="optimise_bond_radius"){
        if (param_value == 0) {
            optimise_bond_radius = false;
        } else {
            optimise_bond_radius = true;
        }
    }
}
