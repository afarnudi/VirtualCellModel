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
                }
            }
        }
    }
    //Tell the class how to inturpret the strings in the config file
    assign_parameters();
    //Check if all the parameters are consistant with the physics/Class enviroment
    consistancy_check();
    //Call the initiliser from the old code
    if (importcoordiantes) {
        import_coordinates(import_file_name);
    } else {
        initialise();
    }
//    initialise(Mesh_file_name);
}
void Chromatin::consistancy_check(){
    if (importcoordiantes) {
        
        if (Num_of_Nodes != 0) {
            string errorMessage = TWARN;
            errorMessage +="The chromatin can be initiated either by \"ImportCoordinates\" or \"GenerateRandomChain\". You cannot select both simultaniously. Consult the template configuration for more information.";
            errorMessage +=TRESET;
            throw std::runtime_error(errorMessage);
        } else {
            ifstream readcoords(import_file_name.c_str());
            if (!readcoords.is_open()) {
                string errorMessage = TWARN;
                errorMessage +="Read Error: Could not read '"+import_file_name+"'";
                errorMessage +=TRESET;
                throw std::runtime_error(errorMessage);
                
            }
        }
    } else if (Num_of_Nodes == 0){
        string errorMessage = TWARN;
        errorMessage +="The chromatin can be initiated either by \"ImportCoordinates\" or \"GenerateRandomChain\". Non was selected. Consult the template configuration for more information.";
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
                    cout<<TWARN<<"Too many arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJsigma\". The first "<<num_of_node_types<<" will be used: ";
                    for (int i=0; i<num_of_node_types; i++) {
                        cout<<TBOLD<<split[i]<<" ";
                        sigma_LJ[i] = stod(split[i]);
                    }
                    cout<<TRESET<<endl;
                } else if (split.size()>num_of_node_types){
                    cout<<TWARN<<"Too few arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJsigma\". Will fill with default values:";
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
                    cout<<TWARN<<"Too many arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJepsilon\". The first "<<num_of_node_types<<" will be used: ";
                    for (int i=0; i<num_of_node_types; i++) {
                        cout<<TBOLD<<split[i]<<" ";
                        epsilon_LJ[i] = stod(split[i]);
                    }
                    cout<<TRESET<<endl;
                    
                } else if (split.size()>num_of_node_types){
                    cout<<TWARN<<"Too few arguments ( Expected "<<TBOLD<<num_of_node_types<<TWARN<<" got "<<split.size()<<TWARN<<") provided for \"LJepsilon\". Will fill with default values:";
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
            if (split[0] != "path/to/my/coordinates.txt") {
                importcoordiantes =true;
                string filepathe = split[0];
                
                import_file_name = split[0];
            }
        } else if (it.first == "GenerateRandomChain") {
            Num_of_Nodes = stoi(split[0]);
        } else if (it.first == "ExportGeneratedCoordinates") {
            if(split[0]=="true"){
                ExportGeneratedCoordinates=true;
            } else if (split[0]=="false"){
                ExportGeneratedCoordinates=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        } else if (it.first == "NodeTypes") {
            num_of_node_types = stoi(split[0]);
        } else if (it.first == "NodeMass") {
            Node_Mass = stod(split[0]);
        } else if (it.first == "NodeRadius") {
            Node_radius = stod(split[0]);
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
            bond_length = stod(split[0]);
        } else if (it.first == "VirtualBondRadius") {
            bond_radius = stod(split[0]);
        } else if (it.first == "OptimiseBondRadius") {
            if(split[0]=="true"){
                optimise_bond_radius=true;
            } else if (split[0]=="false"){
                optimise_bond_radius=false;
            } else {
                string errorMessage = TWARN;
                errorMessage+="I don't understand  \""+split[0]+"\". Use \"true\" or \"false\".";
                errorMessage+= TRESET;
                throw std::runtime_error(errorMessage);
            }
        }
        
    }
}

void Chromatin::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    string resume_file_name;
    string import_file_name;
    ifstream read_config_file(config_file_name.c_str());
    bool resume_flag=false;
    bool import_flag=false;
    GenConst::ChromatinVirtualSites = false;
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
                    cout<<"Resume flag "<<TOFF<<"off"<<TRESET<<".\n";
                } else {
                    resume_flag=true;
                    resume_file_name=split[2];
                    cout<<"Resume flag "<<TON<<"on"<<TRESET<<". Chromatin will resume using the '"<<TFILE<<resume_file_name<<TRESET<<" file.\n";
                }
            } else if(split[0]=="epsilon"){
                if (split.size()<num_of_node_types + 1) {
                    cout<<TWWARN<<"Warning!!!"<<TRESET<<"Too few arguments for the Lenard Jones epsilon interaction of node types.\nNeed "<<num_of_node_types<<" argumens.\n "<<split.size()-1<<" was provided.\nWill resume chromatin interactions with default value, 0.\n";
                    
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
                    cout<<TWWARN<<"Warning!!!"<<TRESET<<"Too few arguments for the Lenard Jones sigma interaction of node types.\nNeed "<<num_of_node_types<<" argumens.\n "<<split.size()-1<<" was provided.\nWill resume chromatin interactions with default value, 1.5x node_radius \n";
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
                    cout<<"Import flag "<<TOFF<<"off"<<TRESET". The Chromatins will be initiated using the config parameters.\n";
                } else {
                    import_flag=true;
                    import_file_name=split[2];
                    cout<<"Import flag "<<TON<<"on"<<TRESET<<". importing coordinatesfrom '"<<TFILE<<import_file_name<<TRESET<<"'.\n";
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
        Shift_position_xyzVector[0]=param_value;
    } else if (param_name=="Shift_in_Y_direction"){
        Shift_position_xyzVector[1]=param_value;
    } else if (param_name=="Shift_in_Z_direction"){
        Shift_position_xyzVector[2]=param_value;
    } else if (param_name=="x_speed"){
        Shift_velocities_xyzVector[0]=param_value;
    } else if (param_name=="y_speed"){
        Shift_velocities_xyzVector[1]=param_value;
    } else if (param_name=="z_speed"){
        Shift_velocities_xyzVector[2]=param_value;
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
