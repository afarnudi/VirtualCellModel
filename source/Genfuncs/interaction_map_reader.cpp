#include "maps.hpp"
#include "Configfile.hpp"
#include "General_constants.h"

using namespace std;

bool is_interaction(string word);

vector<vector<int> > parse_interactiontable_parameters(vector<string> lines){
    INTERindex iindex;
    vector<vector<int> > inter_table;
    int table_size = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs + GenConst::Num_of_Chromatins ;
    inter_table.resize(table_size);
    for (int i=0; i<table_size; i++) {
        inter_table[i].resize(table_size,0);
    }
    
    int row =1;
    //erase the header
    lines.erase(lines.begin());
    
    for (int i=0; i<table_size ; i++) {
        
        string line = lines[i];
        vector<string> split = split_and_check_for_comments(line);
        
        if (split.size()==0) {
            i--;
        }
        if (split.size()>0) {
            if (is_interaction(split[1])) {
                for (int i=0; i<row; i++) {
                    inter_table[row-1][i] = iindex.INTERACTION[split[1+i]];
                }
                row++;
            }
        }
        //I have placed this here for the programme to stop if there is a bigger table in the config file;
        if (row-1 == table_size) {
            break;
        }
    }
    
    return inter_table;
}

bool is_interaction(string word){
    bool verdict = false;
    INTERindex iindex;
    
    if (iindex.INTERACTION.count(word)>0) {
        verdict = true;
    }
    return verdict;
}


void read_interaction_map(vector<vector<int> > &inter_map){
    
    int map_size = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs + GenConst::Num_of_Chromatins ;
    inter_map.resize(map_size);
    for (int i=0; i<map_size; i++) {
        inter_map[i].resize(map_size,0);
    }
    
//    cout<<"GenConst::Interaction_map = "<<GenConst::Interaction_map<<endl;
    if (GenConst::Interaction_map) {
        
        ifstream read_map(GenConst::Interaction_map_file_name.c_str());
        if (read_map.is_open()) {
            if (!GenConst::Testmode) {
                cout<<"'"<<TFILE<<GenConst::Interaction_map_file_name<<TRESET<<"' interaction map file opened "<<TSUCCESS<<"successfully"<<TRESET<<".\n";
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
            cout<<TFAILED<<"Couldn't open"<<TRESET<<" the interaction map file.\n";
            exit(EXIT_FAILURE);
        }
        
    }
    
}
