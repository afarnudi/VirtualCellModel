//
//  ECM_Node_Pair_Identifier.cpp
//  Mem
//
//  Created by Ali Farnudi on 18/08/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "ECM.h"

void ECM::Node_Bond_identifier_2D(void){
    Num_of_Node_Pairs=0;
    
    vector<int> Node_Pairs;
    Node_Pairs.resize(2);
    
//    int ECM_pyramid_Node_1, ECM_pyramid_Node_2, ECM_pyramid_Node_3, ECM_pyramid_Node_4;
//    double temp_double;
//    int Num_of_objects;
    
    int triangle_Node_A, triangle_Node_B, triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Num_of_Triangles;i++)
    {
        triangle_Node_A= Triangle_List[i][0];
        triangle_Node_B= Triangle_List[i][1];
        triangle_Node_C= Triangle_List[i][2];
        
        for(int j=0;j<Node_Bond_list.size();j++)
        {
            if(  ( Node_Bond_list[j][0]==triangle_Node_A &&  Node_Bond_list[j][1]==triangle_Node_B )  ||
                 ( Node_Bond_list[j][0]==triangle_Node_B &&  Node_Bond_list[j][1]==triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( Node_Bond_list[j][0]==triangle_Node_B &&  Node_Bond_list[j][1]==triangle_Node_C )  ||
                 ( Node_Bond_list[j][0]==triangle_Node_C &&  Node_Bond_list[j][1]==triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( Node_Bond_list[j][0]==triangle_Node_A &&  Node_Bond_list[j][1]==triangle_Node_C )  ||
                 ( Node_Bond_list[j][0]==triangle_Node_C &&  Node_Bond_list[j][1]==triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            Node_Pairs[0]=triangle_Node_A;
            Node_Pairs[1]=triangle_Node_B;
            
            Node_Bond_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber2==0)
        {
            Node_Pairs[0]=triangle_Node_B;
            Node_Pairs[1]=triangle_Node_C;
            
            Node_Bond_list.push_back(Node_Pairs);
        }
        
        if(repeatednumber3==0)
        {
            Node_Pairs[0]=triangle_Node_A;
            Node_Pairs[1]=triangle_Node_C;
            
            Node_Bond_list.push_back(Node_Pairs);
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    Num_of_Node_Pairs=int(Node_Bond_list.size());
    cout<<"ECM # of node pairs: "<<Node_Bond_list.size()<<endl;
}






void ECM::Node_Bond_identifier_3D(void){

    vector<int> push;
    push.resize(2);

    //The first Pyramid:
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][1];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][2];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][0];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][1];
    push[1]=Pyramid_Nodes[0][2];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][1];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=Pyramid_Nodes[0][2];
    push[1]=Pyramid_Nodes[0][3];
    Node_Bond_list.push_back(push);

    int pyramid_Node_1, pyramid_Node_2, pyramid_Node_3, pyramid_Node_4;

    bool repeated_pair_1=false;
    bool repeated_pair_2=false;
    bool repeated_pair_3=false;
    bool repeated_pair_4=false;
    bool repeated_pair_5=false;
    bool repeated_pair_6=false;



    for(int i=1;i<Pyramid_Nodes.size();i++)
    {
        pyramid_Node_1=Pyramid_Nodes[i][0];
        pyramid_Node_2=Pyramid_Nodes[i][1];
        pyramid_Node_3=Pyramid_Nodes[i][2];
        pyramid_Node_4=Pyramid_Nodes[i][3];

        for(int j=0;j<Node_Bond_list.size();j++)
        {
            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_2 )  || (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_1=true;
            }

            if(  (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_3 )  || (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_2 )    )
            {
                repeated_pair_2=true;
            }

            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_3 )  || (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_3=true;
            }

            if(  (  Node_Bond_list[j][0] == pyramid_Node_1 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_1 )    )
            {
                repeated_pair_4=true;
            }

            if(  (  Node_Bond_list[j][0] == pyramid_Node_2 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_2 )    )
            {
                repeated_pair_5=true;
            }

            if(  (  Node_Bond_list[j][0] == pyramid_Node_3 &   Node_Bond_list[j][1] == pyramid_Node_4 )  || (  Node_Bond_list[j][0] == pyramid_Node_4 &   Node_Bond_list[j][1] == pyramid_Node_3 )    )
            {
                repeated_pair_6=true;
            }
        }

        if(!repeated_pair_1)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_2;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_2)
        {
            push[0]=pyramid_Node_2;
            push[1]=pyramid_Node_3;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_3)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_3;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_4)
        {
            push[0]=pyramid_Node_1;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_5)
        {
            push[0]=pyramid_Node_2;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }
        if(!repeated_pair_6)
        {
            push[0]=pyramid_Node_3;
            push[1]=pyramid_Node_4;
            Node_Bond_list.push_back(push);
        }

        repeated_pair_1=false;
        repeated_pair_2=false;
        repeated_pair_3=false;
        repeated_pair_4=false;
        repeated_pair_5=false;
        repeated_pair_6=false;

    }

    Num_of_Node_Pairs=int(Node_Bond_list.size());
//    cout<<"# of node pairs: "<<Num_of_Node_Pairs<<endl;
}







void ECM::Node_Bond_identifier_3D_square(void){

    vector<int> push;
    push.resize(2);

    //The first Pyramid:
    push[0]=square_Nodes[0][0];
    push[1]=square_Nodes[0][1];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][0];
    push[1]=square_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][0];
    push[1]=square_Nodes[0][4];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][1];
    push[1]=square_Nodes[0][2];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][1];
    push[1]=square_Nodes[0][5];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][2];
    push[1]=square_Nodes[0][3];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][2];
    push[1]=square_Nodes[0][6];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][3];
    push[1]=square_Nodes[0][7];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][4];
    push[1]=square_Nodes[0][5];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][4];
    push[1]=square_Nodes[0][7];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][5];
    push[1]=square_Nodes[0][6];
    Node_Bond_list.push_back(push);
    push[0]=square_Nodes[0][6];
    push[1]=square_Nodes[0][7];
    Node_Bond_list.push_back(push);

    int square_Node_1, square_Node_2, square_Node_3, square_Node_4;
    int square_Node_5, square_Node_6, square_Node_7, square_Node_8;

    bool repeated_pair_1=false;
    bool repeated_pair_2=false;
    bool repeated_pair_3=false;
    bool repeated_pair_4=false;
    bool repeated_pair_5=false;
    bool repeated_pair_6=false;
    bool repeated_pair_7=false;
    bool repeated_pair_8=false;
    bool repeated_pair_9=false;
    bool repeated_pair_10=false;
    bool repeated_pair_11=false;
    bool repeated_pair_12=false;



    for(int i=1;i<square_Nodes.size();i++)
    {
        square_Node_1=square_Nodes[i][0];
        square_Node_2=square_Nodes[i][1];
        square_Node_3=square_Nodes[i][2];
        square_Node_4=square_Nodes[i][3];
        square_Node_5=square_Nodes[i][4];
        square_Node_6=square_Nodes[i][5];
        square_Node_7=square_Nodes[i][6];
        square_Node_8=square_Nodes[i][7];

        for(int j=0;j<Node_Bond_list.size();j++)
        {
            if(  (  Node_Bond_list[j][0] == square_Node_1 &   Node_Bond_list[j][1] == square_Node_2 )  || (  Node_Bond_list[j][0] == square_Node_2 &   Node_Bond_list[j][1] == square_Node_1 )    )
            {
                repeated_pair_1=true;
            }
            
            if(  (  Node_Bond_list[j][0] == square_Node_1 &   Node_Bond_list[j][1] == square_Node_4 )  || (  Node_Bond_list[j][0] == square_Node_4 &   Node_Bond_list[j][1] == square_Node_1 )    )
            {
                repeated_pair_2=true;
            }
            
            if(  (  Node_Bond_list[j][0] == square_Node_1 &   Node_Bond_list[j][1] == square_Node_5 )  || (  Node_Bond_list[j][0] == square_Node_5 &   Node_Bond_list[j][1] == square_Node_1 )    )
            {
                repeated_pair_3=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_2 &   Node_Bond_list[j][1] == square_Node_3 )  || (  Node_Bond_list[j][0] == square_Node_3 &   Node_Bond_list[j][1] == square_Node_2 )    )
            {
                repeated_pair_4=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_2 &   Node_Bond_list[j][1] == square_Node_6 )  || (  Node_Bond_list[j][0] == square_Node_6 &   Node_Bond_list[j][1] == square_Node_2 )    )
            {
                repeated_pair_5=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_3 &   Node_Bond_list[j][1] == square_Node_4 )  || (  Node_Bond_list[j][0] == square_Node_4 &   Node_Bond_list[j][1] == square_Node_3 )    )
            {
                repeated_pair_6=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_3 &   Node_Bond_list[j][1] == square_Node_7 )  || (  Node_Bond_list[j][0] == square_Node_7 &   Node_Bond_list[j][1] == square_Node_3 )    )
            {
                repeated_pair_7=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_4 &   Node_Bond_list[j][1] == square_Node_8 )  || (  Node_Bond_list[j][0] == square_Node_8 &   Node_Bond_list[j][1] == square_Node_4 )    )
            {
                repeated_pair_8=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_5 &   Node_Bond_list[j][1] == square_Node_6 )  || (  Node_Bond_list[j][0] == square_Node_6 &   Node_Bond_list[j][1] == square_Node_5 )    )
            {
                repeated_pair_9=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_5 &   Node_Bond_list[j][1] == square_Node_8 )  || (  Node_Bond_list[j][0] == square_Node_8 &   Node_Bond_list[j][1] == square_Node_5 )    )
            {
                repeated_pair_10=true;
            }
            
            
            if(  (  Node_Bond_list[j][0] == square_Node_6 &   Node_Bond_list[j][1] == square_Node_7 )  || (  Node_Bond_list[j][0] == square_Node_7 &   Node_Bond_list[j][1] == square_Node_6 )    )
            {
                repeated_pair_11=true;
            }
            
            if(  (  Node_Bond_list[j][0] == square_Node_7 &   Node_Bond_list[j][1] == square_Node_8 )  || (  Node_Bond_list[j][0] == square_Node_8 &   Node_Bond_list[j][1] == square_Node_7 )    )
            {
                repeated_pair_12=true;
            }

            
        }

        if(!repeated_pair_1)
        {
            push[0]=square_Node_1;
            push[1]=square_Node_2;
            Node_Bond_list.push_back(push);
        }
        
        if(!repeated_pair_2)
               {
                   push[0]=square_Node_1;
                   push[1]=square_Node_4;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_3)
               {
                   push[0]=square_Node_1;
                   push[1]=square_Node_5;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_4)
               {
                   push[0]=square_Node_2;
                   push[1]=square_Node_3;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_5)
               {
                   push[0]=square_Node_2;
                   push[1]=square_Node_6;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_6)
               {
                   push[0]=square_Node_3;
                   push[1]=square_Node_4;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_7)
               {
                   push[0]=square_Node_3;
                   push[1]=square_Node_7;
                   Node_Bond_list.push_back(push);
               }
        
        if(!repeated_pair_8)
               {
                   push[0]=square_Node_4;
                   push[1]=square_Node_8;
                   Node_Bond_list.push_back(push);
               }
        
        
        if(!repeated_pair_9)
        {
            push[0]=square_Node_5;
            push[1]=square_Node_6;
            Node_Bond_list.push_back(push);
        }
        
        if(!repeated_pair_10)
        {
            push[0]=square_Node_5;
            push[1]=square_Node_8;
            Node_Bond_list.push_back(push);
        }
        
        if(!repeated_pair_11)
        {
            push[0]=square_Node_6;
            push[1]=square_Node_7;
            Node_Bond_list.push_back(push);
        }
        
        if(!repeated_pair_12)
        {
            push[0]=square_Node_7;
            push[1]=square_Node_8;
            Node_Bond_list.push_back(push);
        }

        
        repeated_pair_1=false;
        repeated_pair_2=false;
        repeated_pair_3=false;
        repeated_pair_4=false;
        repeated_pair_5=false;
        repeated_pair_6=false;
        repeated_pair_7=false;
        repeated_pair_8=false;
        repeated_pair_9=false;
        repeated_pair_10=false;
        repeated_pair_11=false;
        repeated_pair_12=false;

    }

    Num_of_Node_Pairs=int(Node_Bond_list.size());
//    cout<<"# of node pairs: "<<Num_of_Node_Pairs<<endl;
}
