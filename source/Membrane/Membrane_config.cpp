#include <sstream>
#include "Membrane.h"
#include "General_constants.h"
using namespace std;

void Membrane::import_config(string config_file_name){
    map<string, double>::iterator it;
    
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    
    if (read_config_file.is_open()) {
        if (!GenConst::Testmode) {
            cout<<"'"<<TFILE<<config_file_name<<TRESET<<"' file opened "<<TSUCCESS<<"successfully"<<TRESET<<".\n";
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
            
            if (split[0] == comment || (split[0][0]=='/' && split[0][1]=='/')) {
                continue;
            }
            
            param_map[split[0]]=stod(split[1]);
            
            if (split[0]=="Resume") {
                
                if (stoi(split[1])==0) {
                    if (!GenConst::Testmode) {
                        cout<<"Resume flag "<<TRED<<"off"<<TRESET<<". Looking for membrane config parameters.\n";
                    }
                    
                } else {
                    resume=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag "<<TGREEN"on"<<TRESET". Membrane will resume using the '"<<resume_file_name<<"' file.\n";
                }
            } else if (split[0]=="Mesh_file_name") {
                if (split[1]=="1"){
                    mesh_format=1;
                }else if (split[1]=="2"){
                    mesh_format=2;
                }
//                cout<<"Mesh format"<<mesh_format;
                Mesh_file_name=split[2];
                if (!GenConst::Testmode) {
                    cout<<"The '"<<TFILE<<Mesh_file_name<<TRESET<<"' file will be used to initilise the Membrane.\n";
                }
                
            } else {
                set_map_parameter(split[0], param_map[split[0]]);
            }
            
        }//End of while(getline(read_config_file, line))
        
        
    } else {//End of if (read_config_file.is_open())
        cout<<TFAILED<<"Couldn't open"<<TRESET<<" the '"<<TFILE<<config_file_name<<TRESET<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    if(!resume && Mesh_file_name=="non"){
        cout<<"The 'Resume' parameter located in the Membrane config file "<<TFAILED<<"is not set"<<TRESET<<"! Resume should be set to 0 for a membrane initilisation or set to 1 if the membrane is to be imported from a 'resume' file. Please edit the membrane config file and run the programme again.\n\nIn case this is a new run, please provide the meshfile name in the mebrane confi file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume) {
        import(resume_file_name);
    } else {
        it=param_map.find("Mesh_file_name");
        if(it!=param_map.end()){
            initialise(Mesh_file_name);
        }

    }
}

void Membrane::set_map_parameter(string param_name, double param_value){
    
//    map<string, double>::iterator it;
    if (param_name=="Node_Mass") {
        Node_Mass=param_value;
    }else if (param_name=="Node_radius"){
        Node_radius=param_value;
    } else if (param_name=="spring_model"){
        spring_model=param_value;
    } else if (param_name=="Spring_coefficient"){
        Spring_coefficient=param_value;
    } else if (param_name=="Bending_coefficient"){
        Bending_coefficient=param_value;
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
    } else if (param_name=="ext_force"){
        ext_force_model=param_value;
    } else if (param_name=="x_force_constant"){
        kx=param_value;
    } else if (param_name=="y_force_constant"){
        ky=param_value;
    } else if (param_name=="z_force_constant"){
        kz=param_value;
    } else if (param_name=="X_in"){
        X_in=param_value;
    } else if (param_name=="Y_in"){
        Y_in=param_value;
    } else if (param_name=="Z_in"){
        Z_in=param_value;
    } else if (param_name=="position_scale_x"){
        X_scale=param_value;
    } else if (param_name=="position_scale_y"){
        Y_scale=param_value;
    } else if (param_name=="position_scale_z"){
        Z_scale=param_value;
    } else if (param_name=="rescale_factor"){
        rescale_factor=param_value;
    } else if (param_name=="sigma_LJ_12_6"){
        sigma_LJ_12_6=param_value;
    } else if (param_name=="epsilon_LJ_12_6"){
        epsilon_LJ_12_6=param_value;
    } else if(param_name=="Update_radius"){
        
        New_node_radius = param_value;
        
    } else if(param_name=="Begin_update_time_in_Ps"){
        
        Begin_update_time_in_Ps = param_value;
        
    } else if(param_name=="End_update_time_in_Ps"){
        
        End_update_time_in_Ps = param_value;
    } else if(param_name=="Update_nominal_length"){
        
        Update_nominal_length = param_value;
    } else if(param_name=="FENE_min"){
        
        FENE_min = param_value;
    } else if(param_name=="FENE_max"){
        
        FENE_max = param_value;
    } else if(param_name=="FENE_eps"){
        
        FENE_epsilon = param_value;
    } else if(param_name=="FENE_k"){
        
        FENE_k = param_value;
    } else if(param_name=="init_random_rotation"){
        if (int(param_value) == 0) {
            initial_random_rotation_coordinates=false;
        } else{
            initial_random_rotation_coordinates=true;
        }
    }
    
}
