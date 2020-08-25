#include <sstream>
#include "Actin.h"
#include "General_constants.h"
using namespace std;

void Actin::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    string resume_file_name;
    Mesh_file_name="non";
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    int type=1;
    
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
                //                        set_parameter(general_param_map, param_name, param_value);
                //                general_param_map[param_name]=param_value;
                if (stoi(split[1])==0) {
                    cout<<"Resume flag off. Looking for Actin config parameters.\n";
                } else {
                    resume=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag on. Actin will resume using the '"<<resume_file_name<<"' file.\n";
                }
//                break;
            } else if (split[0]=="Mesh_file_name") {
                Mesh_file_name=split[2];
                type=stoi(split[1]);
                cout<<"The '"<<Mesh_file_name<<"' file will be used to initilise the Actin.\n";
                if (type==1) {
                    cout<<"The '"<<Mesh_file_name<<"' will be initialised as a normal type.";
                } else if (type==2){
                    cout<<"The '"<<Mesh_file_name<<"' will be initialised as a 3filaments type.";
                }
//                break;
            } else {
                set_map_parameter(split[0], param_map[split[0]]);
                //                break;
                //                cout<<split[i]<<"\t";
            }
            
        }//End of while(getline(read_config_file, line))
        
        
    } else {//End of if (read_config_file.is_open())
        cout<<"Couldn't open the '"<<config_file_name<<"' file.\n";
        exit(EXIT_FAILURE);
    }
    if(!resume && Mesh_file_name=="non"){
        cout<<"The 'Resume' parameter located in the Membrane config file is not set! Resume should be set to 0 for a membrane initilisation or set to 1 if the membrane is to be imported from a 'resume' file. Please edit the membrane config file and run the programme again.\n\nIn case this is a new run, please provide the meshfile name in the mebrane confi file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume) {
//        import(resume_file_name);
    } else {
        it=param_map.find("Mesh_file_name");
        if(it!=param_map.end()){
            initialise(type);
        }
//        else {
//            cout<<"Resume is off and no meshfile name is provided for initilisation. Please check the membrane config file.\n";
//        }
    }
}

void Actin::set_map_parameter(string param_name, double param_value){
    
//    map<string, double>::iterator it;
    if (param_name=="Node_Mass") {
        Node_Mass=param_value;
    } else if (param_name=="Node_radius"){
        Node_radius=param_value;
    } else if (param_name=="spring_model"){
        spring_model=param_value;
    } else if (param_name=="abp_spring_model"){
        abp_spring_model=param_value;
    } else if (param_name=="MT_spring_model"){
        MT_spring_model=param_value;
    } else if (param_name=="Spring_coefficient"){
        Spring_coefficient=param_value*GenConst::MD_T;
    } else if (param_name=="abp_Spring_coefficient"){
        abp_Spring_coefficient=param_value*GenConst::MD_T;
    } else if (param_name=="MT_Spring_coefficient"){
        MT_Spring_coefficient=param_value*GenConst::MD_T;
    } else if (param_name=="Contractile_force"){
        Contractile_force=param_value;
    } else if (param_name=="Contractile_model"){
        contractile_model=param_value;
    } else if (param_name=="Contractile_hill_co"){
        contractile_hill_co=param_value;
    } else if (param_name=="Contractile_k1"){
        Contractile_k1=param_value;
    } else if (param_name=="Contractile_k2"){
        Contractile_k2=param_value;
    } else if (param_name=="Contractile_rmin_factor"){
        Contractile_rmin=param_value;
    } else if (param_name=="Contractile_rmax_factor"){
        Contractile_rmax=param_value;
    } else if (param_name=="abp_force"){
        abp_force=param_value;
    } else if (param_name=="abp_model"){
        abp_model=param_value;
    } else if (param_name=="abp_hill_co"){
        abp_hill_co=param_value;
    } else if (param_name=="abp_k1"){
        abp_k1=param_value;
    } else if (param_name=="abp_k2"){
        abp_k2=param_value;
    } else if (param_name=="abp_rmin_factor"){
        abp_rmin=param_value;
    } else if (param_name=="abp_rmax_factor"){
        abp_rmax=param_value;
    } else if (param_name=="abp_r0factor"){
        abp_r0factor=param_value;
        
    } else if (param_name=="MT_force"){
            MT_force=param_value;
    } else if (param_name=="MT_model"){
            MT_model=param_value;
    } else if (param_name=="MT_hill_co"){
            MT_hill_co=param_value;
    } else if (param_name=="MT_k1"){
            MT_k1=param_value;
    } else if (param_name=="MT_k2"){
        MT_k2=param_value;
    } else if (param_name=="MT_rmin_factor"){
            MT_rmin=param_value;
    } else if (param_name=="MT_rmax_factor"){
            MT_rmax=param_value;
    } else if (param_name=="MT_r0factor"){
            MT_r0factor=param_value;
        
    } else if (param_name=="actin_r0factor"){
        act_r0factor=param_value;
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
    } else if (param_name=="Kelvin_Damping_Coefficient"){
        Kelvin_Damping_Coefficient=param_value;
    } else if (param_name=="Dashpot_Viscosity"){
        Dashpot_Viscosity=param_value;
    } else if (param_name=="rescale_factor"){
        rescale_factor=param_value;
    }
}
