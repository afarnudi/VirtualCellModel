#include <sstream>
#include "ECM.h"
#include "General_constants.h"
using namespace std;

void ECM::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    string resume_file_name;
    Mesh_file_name="non";
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    int dimension=2;
    
    if (read_config_file.is_open()) {
        cout<<"'"<<TFILE<<config_file_name<<TRESET<<"' file opened "<<TSUCCESS<<"successfully"<<TRESET".\n";
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
                //                        set_parameter(general_param_map, param_name, param_value);
                //                general_param_map[param_name]=param_value;
                if (stoi(split[1])==0) {
                    cout<<"Resume flag "<<TOFF<<"off"<<TRESET<<". Looking for ECM config parameters.\n";
                } else {
                    resume=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag "<<TON<<"on"<<TRESET<<". ECM will resume using the '"<<TFILE<<resume_file_name<<TRESET<<"' file.\n";
                }
            } else if (split[0]=="Mesh_file_name") {
                Mesh_file_name=split[2];
                dimension=stoi(split[1]);
                cout<<"The '"<<TFILE<<Mesh_file_name<<TRESET<<"' file will be used to initilise the ECM.\n";
                if (dimension==2) {
                    cout<<"The '"<<TFILE<<Mesh_file_name<<TRESET<<"' will be initialised as a 2D mesh.";
                } else if (dimension==3){
                    cout<<"The '"<<TFILE<<Mesh_file_name<<TRESET<<"' will be initialised as a 3D mesh.";
                }
            } else {
                set_map_parameter(split[0], param_map[split[0]]);
                //                break;
                //                cout<<split[i]<<"\t";
            }
            
        }//End of while(getline(read_config_file, line))
        
        
    } else {//End of if (read_config_file.is_open())
        cout<<TWWARN<<"Couldn't open the '"<<TFILE<<config_file_name<<TRESET<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    if(!resume && Mesh_file_name=="non"){
        cout<<TWWARN<<"Warning!!!"<<TRESET<<"The 'Resume' parameter located in the Membrane config file is not set! Resume should be set to 0 for a membrane initilisation or set to 1 if the membrane is to be imported from a 'resume' file. Please edit the membrane config file and run the programme again.\n\nIn case this is a new run, please provide the meshfile name in the mebrane confi file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume) {
//        import(resume_file_name);
    } else {
        it=param_map.find("Mesh_file_name");
        if(it!=param_map.end()){
            initialise(dimension);
        }
//        else {
//            cout<<"Resume is off and no meshfile name is provided for initilisation. Please check the membrane config file.\n";
//        }
    }
}

void ECM::set_map_parameter(string param_name, double param_value){
    
//    map<string, double>::iterator it;
    if (param_name=="Node_Mass") {
        Node_Mass=param_value;
    } else if (param_name=="Node_radius"){
        Node_radius=param_value;
    } else if (param_name=="spring_model"){
        spring_model=param_value;
    } else if (param_name=="Spring_coefficient"){
        Spring_coefficient=param_value*GenConst::MD_T;
    } else if (param_name=="receptor_density"){
        receptor_density=param_value;
    } else if (param_name=="receptor_type"){
        receptor_type=param_value;
    }  else if (param_name=="Shift_in_X_direction"){
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
    } else if (param_name=="stiffness_gradient_x"){
        stiffness_gradient_x=param_value;
    } else if (param_name=="stiffness_gradient_y"){
        stiffness_gradient_y=param_value;
    } else if (param_name=="stiffness_gradient_z"){
        stiffness_gradient_z=param_value;
    } else if (param_name=="receptor_gradient_x"){
            receptor_gradient_x=param_value;
    } else if (param_name=="receptor_gradient_y"){
            receptor_gradient_y=param_value;
    } else if (param_name=="receptor_gradient_z"){
            receptor_gradient_z=param_value;
    } else if (param_name=="receptor_center_x"){
                receptor_center_x=param_value;
    } else if (param_name=="receptor_center_y"){
                receptor_center_y=param_value;
    } else if (param_name=="receptor_center_z"){
                receptor_center_z=param_value;
    } else if (param_name=="Kelvin_Damping_Coefficient"){
        Kelvin_Damping_Coefficient=param_value;
    } else if (param_name=="Dashpot_Viscosity"){
        Dashpot_Viscosity=param_value;
    } else if (param_name=="rescale_factor"){
        rescale_factor=param_value;
    } else if (param_name=="sigma_LJ_12_6"){
        sigma_LJ_12_6=param_value;
    } else if (param_name=="epsilon_LJ_12_6"){
        epsilon_LJ_12_6=param_value;
    }
}
