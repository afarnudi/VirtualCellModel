#include <sstream>
#include "Chromatin.h"
#include "Configfile.hpp"

using namespace std;

void Chromatin::import_config(vector<string> configlines){
    //first line is the 'key'
    configlines.erase(configlines.begin());
    //replace the default values with the parameters read from the Config file
    for (int i=0; i<configlines.size(); i++) {
        if(configlines[i].size()!=0){
            vector<string> split = split_and_check_for_comments(configlines[i], "Chromatin: Config reader");
            if (split.size()!=0) {
                map<string, vector<string> >::iterator it;
                it = Params.find(split[0]);
                //Only replace parameters that actually exist in the Chromatin parameters and ignor anythin else;
                if (it != Params.end()) {
                    configlines[i].erase(configlines[i].begin(),configlines[i].begin()+int(split[0].size()) );
                    it->second[0] = configlines[i];
                } else {
                    cout<<TWARN<<"Note: \""<<TFILE<<split[0]<<TWARN<<"\" is not a Chromatin parameter."<<TRESET<<endl;
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
    initialise();
    
//    initialise(Mesh_file_name);
}
void Chromatin::consistancy_check(){
    
    if (Node_radius_stat!="Av") {
        try {
            int test = stoi(Node_radius_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Chromatin config parser: Invalid input for the \"NodeRadius\" (";
            errorMessage+=TFILE;
            errorMessage+=Node_radius_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Av , or just an input integer (example 100)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        Node_radius=stod(Node_radius_stat);
    }
    
    if (BondNominalLength_stat!= "Au" && BondNominalLength_stat!= "Av") {
        try {
            stod(BondNominalLength_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Chromatin config parser: Invalid input for the \"NominalLengthInNm\" (";
            errorMessage+=TFILE;
            errorMessage+=BondNominalLength_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: Au , Av , or just an input value (example 100)";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
    }
    if (BondNominalLength_stat!= "Au") {
        //This will later be used in the OpenMM bond implimentation to define global or per atom parameters. Global will take less memmory.
        spring_model_global_stat=true;
    }
    {
        try {
            stod(SpontaneousBendingAngleInDeg_stat);
        } catch (...) {
            string errorMessage = TWARN;
            errorMessage+="Chromatin config parser: Invalid input for the \"SpontaneousBendingAngleInDeg\" (";
            errorMessage+=TFILE;
            errorMessage+=SpontaneousBendingAngleInDeg_stat;
            errorMessage+=TWARN;
            errorMessage+="). Please try again.\nExample inputs: 180";
            errorMessage+= TRESET;
            throw std::runtime_error(errorMessage);
        }
        SpontaneousBendingAngleInDeg = stod(SpontaneousBendingAngleInDeg_stat);
    }
    
    if (ImportCoordinates) {
        
        if (GenerateRandomChain) {
            string errorMessage = TWARN;
            errorMessage +="The Chromatin can be initiated either by \"ImportCoordinates\" or \"GenerateRandomChain\". You cannot select both simultaniously. Please consult the template configuration for more information.";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        } else {
            ifstream readcoords(import_file_name.c_str());
            if (!readcoords.is_open()) {
                string errorMessage = TWARN;
                errorMessage +="Chromatin config parser: Read Error: Could not read '"+import_file_name+"'";
                errorMessage +=TRESET;
                throw std::runtime_error(errorMessage);
                
            }
        }
    } else if (GenerateRandomChain){
        if (BondNominalLength_stat== "Au" || BondNominalLength_stat== "Av") {
            string errorMessage = TWARN;
            errorMessage +="Chromatin config parser: Cannot generate a random walk path without a step size. Please specify a value for \"NominalLengthInNm\".";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        if (Node_radius_stat== "Av") {
            string errorMessage = TWARN;
            errorMessage +="Chromatin config parser: Cannot generate a random walk. Please specify a value for \"NodeRadius\". This value is used to make sure generated node coordinates do not overlap.";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        }
        
        if (stod(BondNominalLength_stat)<2*stod(Node_radius_stat)) {
            string errorMessage = TWARN;
            errorMessage +="Chromatin config parser: Cannot generate a random walk. Since the \"NodeRadius\" is used to make sure generated node coordinates do not overlap, It cannot be greater than half of the NominalLengthInNm.\n\nNominalLengthInNm = ";
            errorMessage +=TFILE;
            errorMessage +=BondNominalLength_stat;
            errorMessage +=TWARN;
            errorMessage +="\nNodeRadius = ";
            errorMessage +=TFILE;
            errorMessage +=Node_radius_stat;
            errorMessage +=TWARN;
            errorMessage +="\n\nPlease edit the configurations and try again.";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        }
        
    } else {
        string errorMessage = TWARN;
        errorMessage +="The Chromatin can be initiated either by \"ImportCoordinates\" or \"GenerateRandomChain\". Non was selected. Please consult the template configuration for more information.";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    
    if (RestrictGeneratedChainRadius_stat==true && RestrictGeneratedChainRadius < Node_radius ) {
        string errorMessage = TWARN;
        errorMessage +="Chromatin config parser: The Chromatin's node radius cannot be larger than the restricted radius.";
        errorMessage +=TRESET;
        throw std::runtime_error(errorMessage);
    }
    
    
    
    
    
    
    
    epsilon_LJ.resize(num_of_node_types,0);
    sigma_LJ.resize(num_of_node_types,2.5*Node_radius);
    for (auto const& it : Params){
            vector<string> split = split_and_check_for_comments(it.second[0], it.first);
            
        if (it.first == "LJsigma") {
            sigma_LJ.resize(num_of_node_types,2.5*Node_radius);
            if (split.size()==1) {
                if (split[0] != "0") {
                    sigma_LJ[0]= stod(split[0]);
                }
            } else {
                if (split.size()>num_of_node_types) {
                    cout<<TWARN<<"Chromatin config check: LJsigma: Too many arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJsigma\". The first "<<num_of_node_types<<" will be used: ";
                    for (int i=0; i<num_of_node_types; i++) {
                        cout<<TBOLD<<split[i]<<" ";
                        sigma_LJ[i] = stod(split[i]);
                    }
                    cout<<TRESET<<endl;
                } else if (split.size()>num_of_node_types){
                    cout<<TWARN<<"Chromatin config check: LJsigma: Too few arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJsigma\". Will fill with default values:";
                    for (int i=0; i<split.size(); i++) {
                        sigma_LJ[i] = stod(split[i]);
                        cout<<TBOLD<<split[i]<<" ";
                    }
                    for (int i=split.size(); i<num_of_node_types; i++) {
                        cout<<TBOLD<<sigma_LJ[i]<<" ";
                    }
                    cout<<TRESET<<endl;
                }
                
            }
        } else if (it.first == "LJepsilon"){
            epsilon_LJ.resize(num_of_node_types,0);
            if (split.size()==1) {
                if (split[0] != "0") {
                    epsilon_LJ[0]= stod(split[0]);
                }
            } else {
                if (split.size()>num_of_node_types) {
                    cout<<TWARN<<"Chromatin config check: LJepsilon: Too many arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJepsilon\". The first "<<num_of_node_types<<" will be used: ";
                    for (int i=0; i<num_of_node_types; i++) {
                        cout<<TBOLD<<split[i]<<" ";
                        epsilon_LJ[i] = stod(split[i]);
                    }
                    cout<<TRESET<<endl;
                    
                } else if (split.size()>num_of_node_types){
                    cout<<TWARN<<"Chromatin config check: LJepsilon: Too few arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJepsilon\". Will fill with default values:";
                    for (int i=0; i<split.size(); i++) {
                        epsilon_LJ[i] = stod(split[i]);
                        cout<<TBOLD<<split[i]<<" ";
                    }
                    for (int i=split.size(); i<num_of_node_types; i++) {
                        cout<<TBOLD<<epsilon_LJ[i]<<" ";
                    }
                    cout<<TRESET<<endl;
                }
                
            }
        }
        
    }
    
    
}
void Chromatin::assign_parameters(void){
    for (auto const& it : Params){
        vector<string> split = split_and_check_for_comments(it.second[0], "Chroamtin: "+it.first);
        
        if (it.first == "ImportCoordinates") {
            if (split.size()>0) {
                if (split[0]!="path/to/my/coordinates.txt") {
                    import_file_name = split[0];
                    ImportCoordinates =true;
                }
                if (split.size()==2){
                    limit_import = stoi(split[1]);
                }
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: Invalid number of inputs for the \"ImportCoordinates\". Expected 1 or 2 inputs, got (";
                errorMessage+=TFILE;
                errorMessage+=split.size();
                errorMessage+=TWARN;
                errorMessage+="). Please consult the template configuration for more information and try again.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
            
            
        } else if (it.first == "GenerateRandomChain") {
            if (split.size() == 0) {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: Invalid input for the \"GenerateRandomChain\". Expected an integer, got non.";
                errorMessage+="Please try again.";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            } else {
                if (split[0]!="0") {
                    GenerateRandomChain = true;
                    Num_of_Nodes = stoi(split[0]);
                }
                
            }
        } else if (it.first == "RestrictGeneratedChainRadius") {
            if (split[0]!="N") {
                RestrictGeneratedChainRadius_stat=true;
                RestrictGeneratedChainRadius = stod(split[0]);
            }
        } else if (it.first == "NominalLengthInNm") {
            BondNominalLength_stat = split[0];
        } else if (it.first == "SpontaneousBendingAngleInDeg") {
            SpontaneousBendingAngleInDeg_stat = split[0];
        } else if (it.first == "GenerateSelfAvoidingChain") {
            if(split[0]=="true"){
                GenerateSelfAvoidingChain=true;
            } else if (split[0]=="false"){
                GenerateSelfAvoidingChain=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: GenerateSelfAvoidingChain: I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "zebraBonds") {
            if(split[0]=="true"){
                zebraBonds=true;
            } else if (split[0]=="false"){
                zebraBonds=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: zebraBonds: I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "ExportGeneratedCoordinates") {
            if(split[0]=="true"){
                ExportGeneratedCoordinates=true;
            } else if (split[0]=="false"){
                ExportGeneratedCoordinates=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: ExportGeneratedCoordinates: I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "NodeTypes") {
            num_of_node_types = stoi(split[0]);
        } else if (it.first == "NodeMass") {
            Node_Mass = stod(split[0]);
        } else if (it.first == "NodeRadius") {
            Node_radius_stat = split[0];
        } else if (it.first == "SpringModel") {
            if (split[0]=="H") {
                spring_model = potentialModelIndex.Model["Harmonic"];
            } else if (split[0]=="KG") {
                spring_model = potentialModelIndex.Model["KremerGrest"];
            } else if (split[0]=="N") {
                spring_model = potentialModelIndex.Model["None"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: Spring Model: I don't understand the \""+split[0]+"\" Model. Available models: H (Harmonic), N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "AngleSpringModel") {
            if (split[0]=="cos") {
                angle_spring_model = potentialModelIndex.Model["angleCOS"];
            } else if (split[0]=="N") {
                angle_spring_model = potentialModelIndex.Model["None"];
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: Angle Spring Model: I don't understand the \""+split[0]+"\" Model. Available models: cos (cosine potential), N (None).";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "SpringCoeff") {
            Spring_coefficient = stod(split[0]);
        } else if (it.first == "DampingCoeff") {
            Damping_coefficient = stod(split[0]);
        } else if (it.first == "AngleSpringCoeff") {
            Angle_spring_coefficient = stod(split[0]);
        } else if (it.first == "CoordinateTranslateVector") {
            Shift_position_xyzVector.resize(3,0);
            Shift_position_xyzVector[0] = stod(split[0]);
            Shift_position_xyzVector[1] = stod(split[1]);
            Shift_position_xyzVector[2] = stod(split[2]);
        } else if (it.first == "VelocityShiftVector") {
            Shift_velocities_xyzVector.resize(3,0);
            Shift_velocities_xyzVector[0] = stod(split[0]);
            Shift_velocities_xyzVector[1] = stod(split[1]);
            Shift_velocities_xyzVector[2] = stod(split[2]);
        } else if (it.first == "Scale") {
            rescale_factor = stod(split[0]);
        } else if (it.first == "VirtualBondLength") {
            //Need to rewrite this part with the nominal length edition.
//            BondNominalLength_stat = split[0];
//            bond_length = stod(split[0]);
        } else if (it.first == "VirtualBondRadius") {
            bond_radius = stod(split[0]);
        } else if (it.first == "OptimiseBondRadius") {
            if(split[0]=="true"){
                optimise_bond_radius=true;
            } else if (split[0]=="false"){
                optimise_bond_radius=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="Chromatin config parser: OptimiseBondRadius: I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        }
        
    }
}
