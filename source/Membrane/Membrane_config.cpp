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
    
    if (Node_Bond_Nominal_Length_stat!= "Au" && Node_Bond_Nominal_Length_stat!= "Av") {
        try {
            Node_Bond_user_defined_Nominal_Length_in_Nm= stod(Node_Bond_Nominal_Length_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Membrane config parser: Invalid input for the \"NominalLengthInNm\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_Bond_Nominal_Length_stat;
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
            if (extension == "ply" || extension == "msh") {
                mesh_format=extension;
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
                spring_model = potentialModelIndex.Model["Harmonic"];
            } else if (split[0]=="FENE") {
                spring_model = potentialModelIndex.Model["FENE"];
            } else if (split[0]=="N") {
                spring_model = potentialModelIndex.Model["None"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Spring Model: I don't understand the \""+split[0]+"\" Model. Available models: H (Harmonic), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SpringCoeff") {
            Spring_coefficient = stod(split[0]);
        } else if (it.first == "DampingCoeff") {
            Damping_coefficient = stod(split[0]);
        } else if (it.first == "BendingModel") {
            if (split[0]=="N") {
                bending_model = potentialModelIndex.Model["None"];
            } else if (split[0]=="cosine") {
                bending_model = potentialModelIndex.Model["Dihedral"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Membrane config parser: Bending Model: I don't understand the \""+split[0]+"\" Model. Available models: cosine (cosine of dihedral angle), and N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "BendingCoeff") {
            Bending_coefficient = stod(split[0]);
        } else if (it.first == "NominalLengthInNm") {
            Node_Bond_Nominal_Length_stat = split[0];
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
        }
//        else if (it.first == "LJsigma") {
//            sigma_LJ_12_6 = stod(split[0]);
//        }
        else if (it.first == "LJepsilon") {
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
        } else if (it.first == "WantGeometricProps") {
            if(split[0]=="true"){
                WantGeometricProps=true;
            } else if (split[0]=="false"){
                WantGeometricProps=false;
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
