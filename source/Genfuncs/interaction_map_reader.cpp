
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
        vector<string> split = split_and_check_for_comments(line, "Interaction table: Config reader");
        
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


