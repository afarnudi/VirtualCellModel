#include <sstream>
#include "Membrane.h"

using namespace std;

void Membrane::import_config(string config_file_name){
    
    map<string, double>::iterator it;
    string resume_file_name, Mesh_file_name="non";
    ifstream read_config_file(config_file_name.c_str());
    bool resume=false;
    
    if (read_config_file.is_open()) {
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
            
            for (int i=0; i<split.size(); i++) {
                if (split[i] == comment || (split[i][0]=='/' && split[i][1]=='/')) {
                    break;
                }
                
                param_map[split[i]]=stod(split[i+1]);
                if (split[i]=="Resume") {
                    //                        set_parameter(general_param_map, param_name, param_value);
                    //                general_param_map[param_name]=param_value;
                    if (stoi(split[i+1])==0) {
                        cout<<"Resume flag off. Looking for membrane config parameters.\n";
                    } else {
                        resume=true;
                        resume_file_name=split[i+2];
                        cout<<"Resume flag on. Membrane will resume using the '"<<resume_file_name<<"' file.\n";
                    }
                    break;
                } else if (split[i]=="Mesh_file_name") {
                    Mesh_file_name=split[i+1];
                    cout<<"Membrane will initilise via '"<<Mesh_file_name<<"' file.\n";
                    break;
                }
                break;
                //                cout<<split[i]<<"\t";
            }
            
        }//End of while(getline(read_config_file, line))
        
        
    }//End of if (read_config_file.is_open())
    
    if(!resume && Mesh_file_name=="non"){
        cout<<"The 'Resume' parameter located in the Membrane config file is not set! Resume should be set to 0 for a membrane initilisation or set to 1 if the membrane is to be imported from a 'resume' file. Please edit the membrane config file and run the programme again.\n\nIn case this is a new run, please provide the meshfile name in the mebrane confi file.\n";
        exit(EXIT_FAILURE);
    }
    if (resume) {
        import(resume_file_name);
    } else {
        it=param_map.find("Mesh_file_name");
        if(it!=param_map.end()){
            initialise(Mesh_file_name);
        } else {
            cout<<"Resume is off and no meshfile name is provided for initilisation. Please check the membrane config file.\n";
        }
    }
}
