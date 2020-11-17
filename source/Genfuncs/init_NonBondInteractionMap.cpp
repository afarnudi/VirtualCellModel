
#include "Interaction_table.hpp"

using namespace std;

NonBondInteractionMap::NonBondInteractionMap(vector<string> lines){
    inter_config_type["0"]=0;
    inter_config_type_reverse[0]="0";
    inter_type_name[0]="None";
    
    inter_config_type["LJ"]=1;
    inter_config_type_reverse[1]="LJ";
    inter_type_name[1]="Lennard-Jones";
    
    inter_config_type["EV"]=2;
    inter_config_type_reverse[2]="EV";
    inter_type_name[2]="Excluded-Volume";
    
    inter_config_type["3"]=3;
    inter_config_type_reverse[3]="3";
    inter_type_name[3]="Lennard-Jones3";
    
    inter_config_type["4"]=4;
    inter_config_type_reverse[4]="4";
    inter_type_name[4]="Lennard-Jones4";
    
    inter_config_type["5"]=5;
    inter_config_type_reverse[5]="5";
    inter_type_name[5]="Lennard-Jones5";
    
    inter_config_type["LJCS"]=6;
    inter_config_type_reverse[6]="LJCS";
    inter_type_name[6]="Lennard-JonesChromatinSpecial0";

    
    table_size = GenConst::Num_of_Membranes + GenConst::Num_of_Actins + GenConst::Num_of_ECMs + GenConst::Num_of_Chromatins;
    
    inter_table.resize(table_size);
    force_report.resize(table_size);
    force_label.resize(table_size);
    use_max_radius.resize(table_size);
    for (int i=0; i<table_size; i++) {
        inter_table[i].resize(table_size,0);
        force_report[i].resize(table_size,0);
        force_label[i].resize(table_size);
        use_max_radius[i].resize(table_size,0);
    }
    lines.erase(lines.begin());
    
    for (int i=0; i<table_size ; i++) {
        
        string line = lines[i];
        vector<string> split = split_and_check_for_comments(line, "Interaction table: Config reader");
        
        
        
        if (split.size()==0) {
            if (i==lines.size()) {
                string errorMessage = TWARN;
                errorMessage+="Number of interaction table rows does not match the number of classes. Expected ";
                errorMessage+=TFILE;
                errorMessage+=to_string(table_size);
                errorMessage+=TWARN;
                errorMessage+=", got ";
                errorMessage+=TFILE;
                errorMessage+=to_string(rows-1);
                errorMessage+=TWARN;
                errorMessage+=". Please edit the interaction table and try again.";
                errorMessage+=TWARN;
                throw std::runtime_error(errorMessage);
            } else {
                i--;
            }
        }
        
        if (split.size()>0) {
                for (int i=0; i<rows; i++) {
                    if (is_interaction(split[i+1])) {
                        inter_table[rows-1][i] = parse_interacton_type(split[1+i]);
                        
                        if (needs_report(split[1+i])) {
                            force_report[rows-1][i] = true;
                            force_sum++;
                        } else {
                            force_report[rows-1][i] = false;
                        }
                        
                        if (optimise_radius(split[1+i])) {
                            use_max_radius[rows-1][i] = true;
                        } else {
                            use_max_radius[rows-1][i] = false;
                        }
                       
                    } else {
                        string errorMessage = TWARN;
                        errorMessage+="Interaction Table parser: \""+split[i+1]+"\" not defined. Use the template generator for information on how to write the interaction table.";
                        errorMessage+= TRESET;
                        throw std::runtime_error(errorMessage);
                    }
                }
                rows++;
        }
        //I have placed this here for the programme to stop if there is a bigger table in the config file;
        if (rows-1 == table_size) {
            break;
        }
        
    }
    check_force_consistency();
    
};
void NonBondInteractionMap::check_force_consistency(){
    if (force_sum!=0) {
        GenConst::WantForce=true;
    }
    
    
    if (force_sum>31) {
        string errorMessage = TWARN;
        errorMessage+="Force Reporter: cannot report more than 30 forces. This limitation is set by OpenMM. Please edit the interaction table and try again.";
        errorMessage+= TRESET;
        throw std::runtime_error(errorMessage);
    }
}

int NonBondInteractionMap::setForceGroup(int row, int col){
    ForceGroupCount++;
    forceGroupLabel[ForceGroupCount]=force_label[row][col]+"_"+inter_config_type_reverse[inter_table[row][col]];
    forceGroupIndex[ForceGroupCount]=ForceGroupCount;
    return forceGroupIndex[ForceGroupCount];
}

int NonBondInteractionMap::setForceGroup(int row, int col, int chromotype){
    ForceGroupCount++;
    forceGroupLabel[ForceGroupCount]=force_label[row][col]+"_"+inter_config_type_reverse[inter_table[row][col]]+to_string(chromotype);
    forceGroupIndex[ForceGroupCount]=ForceGroupCount;
    return forceGroupIndex[ForceGroupCount];
}


int NonBondInteractionMap::get_ForceGroup(int index){
    int forcegroup = 1 << forceGroupIndex[index];
    return forcegroup;
}

bool NonBondInteractionMap::is_interaction(string word){
    bool verdict = false;
    string ups_digits;
    for (auto &a: word){
        if (isupper(a) || isdigit(a)) {
            ups_digits+=a;
        }
    }
    if (inter_config_type.count(ups_digits)>0) {
        verdict = true;
    }
    return verdict;
}

bool NonBondInteractionMap::needs_report(string word){
    bool verdict = false;
    string lows;
    for (auto &a: word){
        if (a == 'r') {
            verdict = true;
        }
    }
    return verdict;
}

bool NonBondInteractionMap::optimise_radius(string word){
    bool verdict = false;
    string lows;
    for (auto &a: word){
        if (a == 'm') {
            verdict = true;
        }
    }
    return verdict;
}


int NonBondInteractionMap::parse_interacton_type(string word){
    string interaction;
    for (auto &a: word){
        if (isupper(a) || isdigit(a)) {
            interaction+=a;
        }
    }
    return inter_config_type[interaction];
}


bool NonBondInteractionMap::get_radius_optimisation_status(int row, int col){
    return use_max_radius[row][col];
}
