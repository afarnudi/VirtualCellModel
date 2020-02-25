#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"


const int EndOfList=-1;

//std::vector< std::pair< int, int > > exclusion_list_generator(Bonds*      bonds,
//                                                              std::string label_1,
//                                                              std::string label_2){
//
//    std::vector< std::pair< int, int > > exclude_bonds;
//    //if shared node --> acin or membrane label
//    for (int i_b=0; bonds[i_b].type != EndOfList; ++i_b) {
//
//        if ( (bonds[i_b].class_label == label_1 + label_2)  || (bonds[i_b].class_label == label_2 + label_1) ) {
//            std::pair< int, int > temp;
//            temp.first=bonds[i_b].atoms[0];
//            temp.second=bonds[i_b].atoms[1];
//            exclude_bonds.push_back(temp);
//        }
//
//
//    }
//    return exclude_bonds;
//}


std::vector< std::pair< int, int > > exclusion_list_generator(Bonds*      bonds){
    
    std::vector< std::pair< int, int > > exclude_bonds;
    //if shared node --> acin or membrane label
    for (int i_b=0; bonds[i_b].type != EndOfList; ++i_b) {
            std::pair< int, int > temp;
            temp.first =bonds[i_b].atoms[0];
            temp.second=bonds[i_b].atoms[1];
            exclude_bonds.push_back(temp);
    }
    return exclude_bonds;
}


void add_exclusion(OpenMM::CustomNonbondedForce* custom_bond,
                   std::vector< std::pair< int, int > > exclude_list)
{
    for (int i=0; i<exclude_list.size(); ++i)
    {
        custom_bond->addExclusion(exclude_list[i].first, exclude_list[i].second);
    }
}
