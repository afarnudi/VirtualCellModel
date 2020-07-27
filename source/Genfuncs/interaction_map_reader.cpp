#include "maps.hpp"
#include "General_constants.h"

using namespace std;

void read_interaction_map(vector<vector<int> > &inter_map){
    
    int map_size = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs + GenConst::Num_of_Chromatins + GenConst::Num_of_pointparticles;
    inter_map.resize(map_size);
    for (int i=0; i<map_size; i++) {
        inter_map[i].resize(map_size,0);
    }
    
//    cout<<"GenConst::Interaction_map = "<<GenConst::Interaction_map<<endl;
    if (GenConst::Interaction_map) {
        
        ifstream read_map(GenConst::Interaction_map_file_name.c_str());
        if (read_map.is_open()) {
            if (!GenConst::Testmode) {
                cout<<"\""<<GenConst::Interaction_map_file_name<<"\" interaction map file opened successfully.\n";
            }
            
            string line;
            int line_num=0;
            string comment="//";
            //        char delimiter=' ';
            while(getline(read_map, line)){
                
                if(line.empty()){
                    continue;
                }
                
                istringstream iss(line);
                vector<string> split(istream_iterator<string>{iss}, istream_iterator<string>());
                
                if (line_num < map_size) {
//                    cout<<line_num<<": ";
                    for (int i=0; i<split.size(); i++) {
                        
//                        cout<<split[i]<<"\t";
                        if (split[i] == comment || (split[i][0]=='/' && split[i][1]=='/')) {
                            break;
                        }
//                        cout<<split[i]<<"\t";
                        if (i<line_num+1) {
//                            cout<<split[i+1]<<endl;
                            inter_map[line_num][i]=stoi(split[i+1]);
                            inter_map[i][line_num]=stoi(split[i+1]);
                            if (stoi(split[i+1])==2) {
                                GenConst::Excluded_volume_interaction=true;
                            }
                        }
                        
                    } // End of for (int i=0; i<split.size(); i++) {
                    line_num++;
//                    cout<<"\n";
                }
            } //End of while(getline(read_config_file, line)){
        } else {
            cout<<"Couldn't open the interaction map file.\n";
            exit(EXIT_FAILURE);
        }
        
    }
    
}
