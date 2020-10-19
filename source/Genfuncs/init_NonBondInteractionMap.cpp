#include "Interaction_table.hpp"

using namespace std;
NonBondInteractionMap::NonBondInteractionMap(vector<string> lines){
    inter_config_type["0"]=0;
    inter_type_name[0]="None";
    inter_config_type["LJ"]=1;
    inter_type_name[1]="Lennard-Jones";
    inter_config_type["EV"]=2;
    inter_type_name[2]="Excluded-Volume";
    inter_config_type["3"]=3;
    inter_type_name[3]="Lennard-Jones3";
    inter_config_type["4"]=4;
    inter_type_name[4]="Lennard-Jones4";
    inter_config_type["5"]=5;
    inter_type_name[5]="Lennard-Jones5";
    inter_config_type["LJCS"]=6;
    inter_type_name[6]="Lennard-JonesChromatinSpecial0";
    
    inter_config_type["LJR"]=1001;
    inter_config_type["EVR"]=1002;
    inter_config_type["LJCSR"]=1006;
    
    table_size = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs + GenConst::Num_of_Chromatins;
    
    inter_table.resize(table_size);
    force_report.resize(table_size);
    for (int i=0; i<table_size; i++) {
        inter_table[i].resize(table_size,0);
        force_report[i].resize(table_size,0);
    }
    lines.erase(lines.begin());
    for (int i=0; i<table_size ; i++) {
        
        string line = lines[i];
        vector<string> split = split_and_check_for_comments(line, "Interaction table: Config reader");
        
        if (split.size()==0) {
            i--;
        }
        if (split.size()>0) {
            
                for (int i=0; i<row; i++) {
                    if (is_interaction(split[i+1])) {
                        inter_table[row-1][i] = inter_config_type[split[1+i]]%1000;
                        if (inter_config_type[split[1+i]]>1000) {
                            force_report[row-1][i] = true;
                            force_sum++;
                        } else {
                            force_report[row-1][i] = false;
                        }
                    } else {
                        string errorMessage = TWARN;
                        errorMessage+="Interaction Table parser: \""+split[i+1]+"\" not defined. Use the template generator for information on how to write the interaction table.";
                        errorMessage+= TRESET;
                        throw std::runtime_error(errorMessage);
                    }
                }
                row++;
        }
        //I have placed this here for the programme to stop if there is a bigger table in the config file;
        if (row-1 == table_size) {
            break;
        }
        check_force_consistency();
    }
    
    
}
void NonBondInteractionMap::check_force_consistency(){
    if (force_sum>32) {
        string errorMessage = TWARN;
        errorMessage+="Force Reporter: cannot report more than 31 forces. This limitation is set by OpenMM. Please edit the interaction table and try again.";
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
}
