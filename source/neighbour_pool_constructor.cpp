//
//  neighbour_pool_constructor.cpp
//  Mem
//
//  Created by Ali Farnudi on 21/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "interaction.hpp"

//void neighbour_pool_constructor (vector<int> &neighbour_pool, Membrane membrane){
//    
//    for (int i=0; i<neighbour_pool.size()-1; i++) {
//        int temp_node=neighbour_pool[i];
//        for (int j=i+1; j<neighbour_pool.size(); j++) {
//            if (neighbour_pool[j]==temp_node) {
//                neighbour_pool.erase(neighbour_pool.begin()+j);
//                j--;
//            }
//        }
//    }
//    
//    int pool_size=int(neighbour_pool.size());
////    cout<<"pool_size="<<pool_size<<endl;
//    
//    for (int i=0; i<pool_size; i++) {
//        int temp_node=neighbour_pool[i];
////        cout<<"neighbour_pool["<<i<<"]="<<neighbour_pool[i]<<endl;
//        for (int k=0; k<membrane.Node_neighbour_list[temp_node].size(); k++) {
//            neighbour_pool.push_back(membrane.Node_neighbour_list[temp_node][k]);
////            cout<<"membrane.node_neighbour_list[temp_node]["<<k<<"]="<<membrane.node_neighbour_list[temp_node][k]<<endl;
//        }
//    }
//    
//    for (int i=0; i<neighbour_pool.size()-1; i++) {
//        int temp_node=neighbour_pool[i];
//        for (int j=i+1; j<neighbour_pool.size(); j++) {
//            if (neighbour_pool[j]==temp_node) {
//                neighbour_pool.erase(neighbour_pool.begin()+j);
//                j--;
//            }
//        }
//    }
//    cout<<"neighbour_pool.size()="<<neighbour_pool.size()<<endl;
//}
