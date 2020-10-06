#include <sstream>
#include "Membrane.h"
#include "Configfile.hpp"
#include "General_constants.h"
using namespace std;

void Membrane::import_config(vector<string> configlines){
    
    //first line is the 'key'
    configlines.erase(configlines.begin());
    
    //replace the default values with the parameters read from the Config file
    for (int i=0; i<configlines.size(); i++) {
        if(configlines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(configlines[i], "Membrane: Config reader");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = Params.find(split[0]);
                //Only replace parameters that actually exist in the Membrane parameters and ignor anythin else;
                if (it != Params.end()) {
                    configlines[i].erase(configlines[i].begin(),configlines[i].begin()+int(split[0].size()) );
                    it->second[0] = configlines[i];
                } else {
                    cout<<TWARN<<"Note: \""<<TFILE<<split[0]<<TWARN<<"\" is not a Membrane parameter."<<TRESET<<endl;
                    cout<<"If you wish to edit the configfile, exit. If not, press any key to continue."<<endl;
                    getchar();
                    cout<<TRESET;
                }
            }
        }
    }
    //Tell the class how to inturpret the strings in the config file
    assign_parameters();
    //Check if all the parameters are consistant with the physics/Class enviroment
    consistancy_check();
    //Call the initiliser from the old code
    initialise(Mesh_file_name);
}

void Membrane::consistancy_check(){
    ifstream readmesh(Mesh_file_name.c_str());
    if (!readmesh.is_open()) {
        string errorMessage = TWARN;
        errorMessage +="Read Error: Could not read '"+Mesh_file_name+"'";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    
    if (Node_Bond_distances_stat!= "Au" && Node_Bond_distances_stat!= "Av") {
        try {
            Node_Bond_Nominal_Length_in_Nm= stod(Node_Bond_distances_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \"NominalLengthInNm\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_Bond_distances_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 2.678)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }
    
    if (Triangle_pair_angle_stat!= "Au" && Triangle_pair_angle_stat!= "Av") {
        try {
            Triangle_pair_Nominal_angle_in_degrees= stod(Triangle_pair_angle_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \"SpontaneousTriangleBendingAngleInDegrees\" (";
            errorMessage+=TFILE;
            errorMessage+=Triangle_pair_angle_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 180)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }
    if (Node_radius_stat!= "Au" && Node_radius_stat!= "Av") {
        try {
            double test = stod(Node_radius_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \"NodeRadius\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_radius_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 100)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }
}

void Membrane::assign_parameters(void){
    for (auto const& it : Params){
        vector<string> split = split_and_check_for_comments(it.second[0], "Membrane: "+it.first);
        if (split.size()==0) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \""+it.first+"\". Value in configuration file";
            errorMessage+=TFILE;
            errorMessage+=it.second[0];
            errorMessage+=TWARN;
            errorMessage+=". Please check the template and try again.";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        if (it.first == "MeshFile") {
            string extension = split[0];
            extension.erase(extension.begin(),extension.begin()+extension.find('.')+1);
            if (extension == "ply") {
                mesh_format=2;
            } else if (extension == "msh"){
                mesh_format=1;
            } else{
                string errorMessage = TWARN;
                errorMessage+="I don't understand the \""+extension+"\" format. Please use the Blender (ply) or Gmesh 2 (msh).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            Mesh_file_name = split[0];
        } else if (it.first == "NodeMass") {
            try {
                Node_Mass = stod(split[0]);
            } catch (...) {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Invalid input for the \"NodeMass\", expected a number but got \"";
                errorMessage+=TFILE;
                errorMessage+=split[0];
                errorMessage+=TWARN;
                errorMessage+="\". Please try again.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "NodeRadius") {
            Node_radius_stat = split[0];
        } else if (it.first == "SpringModel") {
            if (split[0]=="H") {
                spring_model = 2;
            } else {
                string errorMessage = TWARN;
                errorMessage+="I don't understand the \""+split[0]+"\" Model. Available models: H (Harmonic).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SpringCoeff") {
            Spring_coefficient = stod(split[0]);
        } else if (it.first == "DampingCoeff") {
            Damping_coefficient = stod(split[0]);
        } else if (it.first == "BendingCoeff") {
            Bending_coefficient = stod(split[0]);
        } else if (it.first == "NominalLengthInNm") {
            Node_Bond_distances_stat = split[0];
        } else if (it.first == "SpontaneousTriangleBendingAngleInDegrees") {
            Triangle_pair_angle_stat = split[0];
        } else if (it.first == "ExtForceModel") {
            ext_force_model = stoi(split[0]);
        } else if (it.first == "XYZinMembrane") {
            X_in = stod(split[0]);
            Y_in = stod(split[1]);
            Z_in = stod(split[2]);
        } else if (it.first == "XYZscale") {
            X_scale = stod(split[0]);
            Y_scale = stod(split[1]);
            Z_scale = stod(split[2]);
        } else if (it.first == "Scale") {
                rescale_factor = stod(split[0]);
        } else if (it.first == "LJsigma") {
            sigma_LJ_12_6 = stod(split[0]);
        } else if (it.first == "LJepsilon") {
            epsilon_LJ_12_6 = stod(split[0]);
        } else if (it.first == "UpdateRadius") {
            New_Radius = stod(split[0]);
        } else if (it.first == "UpdateBeginTimeInPs") {
            Begin_update_time_in_Ps = stod(split[0]);
        } else if (it.first == "UpdateEndTimeInPs") {
            End_update_time_in_Ps = stod(split[0]);
        } else if (it.first == "InitRandomRotation") {
            if(split[0]=="true"){
                initial_random_rotation_coordinates=true;
            } else if (split[0]=="false"){
                initial_random_rotation_coordinates=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "VelocityShiftVector") {
            Shift_velocities_xyzVector.resize(3,0);
            Shift_velocities_xyzVector[0] = stod(split[0]);
            Shift_velocities_xyzVector[1] = stod(split[1]);
            Shift_velocities_xyzVector[2] = stod(split[2]);
        } else if (it.first == "CoordinateTranslateVector") {
            Shift_position_xyzVector.resize(3,0);
            Shift_position_xyzVector[0] = stod(split[0]);
            Shift_position_xyzVector[1] = stod(split[1]);
            Shift_position_xyzVector[2] = stod(split[2]);
        }
        
    }
}

void Membrane::import_config(string config_file_name){
    map<string, double>::iterator it;
    
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    
    Shift_position_xyzVector.resize(3,0);
    Shift_velocities_xyzVector.resize(3,0);
    
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
                        cout<<"Resume flag "<<TOFF<<"off"<<TRESET<<". Looking for membrane config parameters.\n";
                    }
                    
                } else {
                    resume=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag "<<TON"on"<<TRESET". Membrane will resume using the '"<<resume_file_name<<"' file.\n";
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
                
            } else if (split[0]=="Shift_positions_xyzVector") {
                Shift_position_xyzVector[0] = stod(split[1]);
                Shift_position_xyzVector[1] = stod(split[2]);
                Shift_position_xyzVector[2] = stod(split[3]);
                
            } else if (split[0]=="Shift_velocities_xyzVector") {
                Shift_velocities_xyzVector[0] = stod(split[1]);
                Shift_velocities_xyzVector[1] = stod(split[2]);
                Shift_velocities_xyzVector[2] = stod(split[3]);
                
                if (!GenConst::Testmode) {
                    cout<<"Shift_velocities_xyzVector "<<Shift_velocities_xyzVector[0]<<" "<<Shift_velocities_xyzVector[1]<<" "<<Shift_velocities_xyzVector[2]<<endl;
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
        Node_radius_stat=param_value;
    } else if (param_name=="spring_model"){
        spring_model=param_value;
    } else if (param_name=="Spring_coefficient"){
        Spring_coefficient=param_value;
    } else if (param_name=="Bending_coefficient"){
        Bending_coefficient=param_value;
    } else if (param_name=="Damping_coefficient"){
        Damping_coefficient=param_value;
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
        New_Radius = param_value;
    } else if(param_name=="Begin_update_time_in_Ps"){
        Begin_update_time_in_Ps = param_value;
    } else if(param_name=="End_update_time_in_Ps"){
        End_update_time_in_Ps = param_value;
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
