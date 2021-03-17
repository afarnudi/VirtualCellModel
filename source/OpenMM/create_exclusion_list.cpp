#include "Membrane.h"
#include "General_functions.hpp"
#include "OpenMM_structs.h"
#include "OpenMM_funcs.hpp"


const int EndOfList=-1;


void creatBondExclusion(Bonds*                                 bonds,
                        NonBondInteractionMap                 &interaction_map,
                        vector<OpenMM::CustomNonbondedForce*> &LJ_12_6_interactions,
                        vector<OpenMM::CustomNonbondedForce*> &ExcludedVolumes,
                        vector<OpenMM::CustomNonbondedForce*> &WCAs
                        ){
    std::vector< std::pair< int, int > > excludedbonds = exclusion_list_generator(bonds);
    
    if (WCAs.size()==1) {
        WCAs[0]->createExclusionsFromBonds(excludedbonds, 1);
//        cout<<WCAs[0]->getNumExclusions()<<endl;
//        exit(0);
    }
    
}


std::vector< std::pair< int, int > > exclusion_list_generator(Bonds*      bonds){
    
    std::vector< std::pair< int, int > > exclude_bonds;
    //if shared node --> acin or membrane label
    for (int i_b=0; bonds[i_b].type != EndOfList; ++i_b) {
            std::pair< int, int > temp;
            temp.first =bonds[i_b].atoms[0];
            temp.second=bonds[i_b].atoms[1];

        int j=0;
        bool flag = true;
        
        while (j<exclude_bonds.size())
        {
            if( (exclude_bonds[j] == temp)  )
            {
                flag = false;
                j=i_b;
            }
            j++;
        }
        
        if(flag)
        {
//            cout<<temp.first<<" "<<temp.second<<endl;
//            if (i_b%20==0) {
//                int a;
//                cin>>a;
//            }
            exclude_bonds.push_back(temp);
        }
        
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
