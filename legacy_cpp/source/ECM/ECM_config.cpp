#include <sstream>
#include "ECM.h"
#include "General_constants.h"
#include "Configfile.hpp"
using namespace std;


void ECM::import_config(vector<string> configlines){
    //first line is the 'key'
    configlines.erase(configlines.begin());
    
    //replace the default values with the parameters read from the Config file
    for (int i=0; i<configlines.size(); i++) {
        if(configlines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(configlines[i], "ECM: Config reader");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = Params.find(split[0]);
                //Only replace parameters that actually exist in the Membrane parameters and ignor anythin else;
                if (it != Params.end()) {
                    configlines[i].erase(configlines[i].begin(),configlines[i].begin()+int(split[0].size()) );
                    it->second[0] = configlines[i];
                } else {
                    cout<<TWARN<<"Note: \""<<TFILE<<split[0]<<TWARN<<"\" is not a ECM parameter."<<TRESET<<endl;
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
  // initialise(Mesh_file_name);
}


void ECM::consistancy_check(){
    ifstream readmesh(Mesh_file_name.c_str());
    if (!readmesh.is_open()) {
        string errorMessage = TWARN;
        errorMessage +="Read Error: Could not read/locate mesh file '"+Mesh_file_name+"'";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }

    if (Node_Bond_Nominal_Length_stat!= "Au" && Node_Bond_Nominal_Length_stat!= "Av") {
        try {
            Node_Bond_user_defined_Nominal_Length_in_Nm= stod(Node_Bond_Nominal_Length_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="ECM config parser: Invalid input for the \"NominalLengthInNm\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_Bond_Nominal_Length_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 2.678)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }

    if (Node_radius_stat!= "Au" && Node_radius_stat!= "Av") {
        try {
            double test = stod(Node_radius_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="ECM config parser: Invalid input for the \"NodeRadius\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_radius_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 100)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }

    if (Spring_coefficient==0) {
        spring_model=potentialModelIndex.Model["None"];
    }

}

void ECM::assign_parameters(void){
    for (auto const& it : Params){
        vector<string> split = split_and_check_for_comments(it.second[0], "ECM: "+it.first);
        if (split.size()==0) {
            string errorMessage = TWARN;
            errorMessage+="ECM config parser: Invalid input for the \""+it.first+"\". Value in configuration file";
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
            if (extension == "ply" || extension == "msh") {
                mesh_format=extension;
            } else{
                string errorMessage = TWARN;
                errorMessage+="I don't understand the \""+extension+"\" format. Please use the Blender (ply) or Gmesh 2 (msh).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            Mesh_file_name = split[0];
            } else if (it.first == "MeshType") {
            MeshType = split[0];
                //3d or 2d meshtype for ecm
        } else if (it.first == "NodeMass") {
            try {
                Node_Mass = stod(split[0]);
            } catch (...) {
                string errorMessage = TWARN;
                errorMessage+="ECM config parser: Invalid input for the \"NodeMass\", expected a number but got \"";
                errorMessage+=TFILE;
                errorMessage+=split[0];
                errorMessage+=TWARN;
                errorMessage+="\". Please try again.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "NodeRadius") {
            Node_radius_stat = split[0];
        } else if (it.first == "NominalLengthInNm") {
            Node_Bond_Nominal_Length_stat = split[0];
        } else if (it.first == "SpringModel") {
            if (split[0]=="H") {
                spring_model = potentialModelIndex.Model["Harmonic"];
            } else if (split[0]=="N") {
                spring_model = potentialModelIndex.Model["None"];
            } else if (split[0]=="kelvin") {
                spring_model = potentialModelIndex.Model["Kelvin-Voigt"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="ECM config parser: Spring Model: I don't understand the \""+split[0]+"\" Model. Available models: H (Harmonic), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SpringCoeff") {
            Spring_coefficient = stod(split[0]);
        } else if (it.first == "DampingCoeff") {
            Kelvin_Damping_Coefficient = stod(split[0]);
        } else if (it.first == "ExtForceModel") {
            ext_force_model = stoi(split[0]);
        } else if (it.first == "Scale") {
                rescale_factor = stod(split[0]);
        } else if (it.first == "stiffness_gradient_x") {
                stiffness_gradient_x = stoi(split[0]);
        } else if (it.first == "stiffness_gradient_y") {
                stiffness_gradient_y = stoi(split[0]);
        } else if (it.first == "stiffness_gradient_z") {
                stiffness_gradient_z = stoi(split[0]);
        } else if (it.first == "receptor_type") {
                receptor_type = stoi(split[0]);
        } else if (it.first == "receptor_density") {
                receptor_density = stoi(split[0]);
        } else if (it.first == "receptor_gradient_x") {
                receptor_gradient_x = stoi(split[0]);
        } else if (it.first == "receptor_gradient_y") {
                receptor_gradient_y = stoi(split[0]);
        } else if (it.first == "receptor_gradient_z") {
                receptor_gradient_z = stoi(split[0]);
        } else if (it.first == "receptor_center_x") {
                receptor_center_x = stoi(split[0]);
        } else if (it.first == "receptor_center_y") {
                receptor_center_y = stoi(split[0]);
        } else if (it.first == "receptor_center_z") {
                receptor_center_z = stoi(split[0]);
        }
//        else if (it.first == "LJsigma") {
//            sigma_LJ_12_6 = stod(split[0]);
//        }
        else if (it.first == "LJepsilon") {
            epsilon_LJ_12_6 = stod(split[0]);
        }
        else if (it.first == "ExtForceModel") {
             ext_force_model = stoi(split[0]);
        
         } else if (it.first == "ExtForceRigidity") {
             ext_force_rigidity.resize(3,0);
             ext_force_rigidity[0] = stoi(split[0]);
             ext_force_rigidity[1] = stoi(split[1]);
             ext_force_rigidity[2] = stoi(split[2]);
        
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
