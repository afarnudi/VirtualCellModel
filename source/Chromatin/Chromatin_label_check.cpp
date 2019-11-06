#include "Chromatin.h"
#include <algorithm>

using std::cout;
using std::endl;

void Chromatin::pdb_label_check(void){
    
    while (label.length()>3) {
        label.pop_back();
    }
    while (label.length()<3) {
        label += "0";
    }
    
    if (index>=100) {
        cout<<"Trouble with PDB labeling:\nPDB name: class_label + calss_index + class_node_type. \nToo many classes (100>) to make a PDB name with the correct format.\n";
        std::exit(EXIT_FAILURE);
    } else if (index>=10){
        if(num_of_node_types>=10){
            cout<<"Trouble with PDB labeling:\nPDB name: class_label + calss_index + class_node_type. \nToo many node types (# of classes >10 and # of node types >10) to make a PDB name with the correct format.\n";
            std::exit(EXIT_FAILURE);
        } else {
            label.pop_back();
            label.pop_back();
            label += std::to_string(index);
        }
    } else {
        if(num_of_node_types>=10){
            label.pop_back();
            label.pop_back();
            label += std::to_string(index);
        } else {
            label.pop_back();
            label += std::to_string(index);
        }
    }
    
    cout<<"label : "<<label<<endl;
}
