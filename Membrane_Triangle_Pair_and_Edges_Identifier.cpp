
#include "Membrane.h"


void Membrane::Membrane_Triangle_Pair_and_Edges_Identifier()
{
    
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C, temp_triangle_node_D, neighbor=0, neighbor_indicator;
    int triangle_pairs=0;
    int temp[2*Membrane_num_of_Triangle_Pairs][4];
    int temp2[2*Membrane_num_of_Triangle_Pairs][4];

	Membrane_Edges.resize(Membrane_num_of_Triangle_Pairs);
	Membrane_Triangle_Pair_Nodes.resize(Membrane_num_of_Triangle_Pairs);
	membrane_triangle_pair_list.resize(Membrane_triangle_list.size());
	for (int i=0; i<Membrane_num_of_Triangle_Pairs; i++)
	{
		Membrane_Edges[i].resize(2);
		Membrane_Triangle_Pair_Nodes[i].resize(4);
	}

	
    
    for(int i=0 ;i<Membrane_triangle_list.size();i++)
    {
        temp_triangle_node_A=Membrane_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        temp_triangle_node_B=Membrane_triangle_list[i][1];
        temp_triangle_node_C=Membrane_triangle_list[i][2];
        neighbor_indicator=0;
		
        //        Indicates the existence of Node neighbor for a node pair (other than the membrane of the triangle):
        //        neighbor_indicator=0 No Node Pairs; neighbor_indicator=1, for temp_triangle_node_A-temp_triangle_node_B; neighbor_indicator=2, for temp_triangle_node_B-temp_triangle_node_C; And neighbor_indicator=3 for temp_triangle_node_C-temp_triangle_node_A.
        for(int j=0;j<Membrane_triangle_list.size();j++)
        {
			
            //************************** finding neighbors **************************
            // neibours of temp_triangle_node_A-temp_triangle_node_B:
            if ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &&  Membrane_triangle_list[j][1]==temp_triangle_node_B  && Membrane_triangle_list[j][2]!=temp_triangle_node_C ){
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=1;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &&  Membrane_triangle_list[j][1]==temp_triangle_node_A  && Membrane_triangle_list[j][2]!=temp_triangle_node_C )
            {
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=1;
            }
			
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=2;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B )
            {
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=2;
            }
			
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=3;
                
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbor=Membrane_triangle_list[j][2];
                neighbor_indicator=3;
                
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbor=Membrane_triangle_list[j][1];
                neighbor_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbor=Membrane_triangle_list[j][0];
                neighbor_indicator=3;
            }
			
			
            
            if(neighbor_indicator!=0)  //  to speed up  the programme we first check if we have found a neighbor or not
            {
                // note that temp_triangle_node_A-temp_triangle_node_B-temp_triangle_node_C-neighbor  are 4 point of two triangle wich will interact
                membrane_triangle_pair_list[i].push_back(j);
                
                if(neighbor_indicator==1)
                {
                    temp[triangle_pairs][0] =temp_triangle_node_A;
                    temp[triangle_pairs][1] =temp_triangle_node_B;
                    temp[triangle_pairs][2] =temp_triangle_node_C;
                    temp[triangle_pairs][3] =neighbor;
                    
                }
                if(neighbor_indicator==2)
                {
                    
                    temp[triangle_pairs][0] =temp_triangle_node_B;
                    temp[triangle_pairs][1] =temp_triangle_node_C;
                    temp[triangle_pairs][2] =temp_triangle_node_A;
                    temp[triangle_pairs][3] =neighbor;
                }
                if(neighbor_indicator==3)
                {
                    
                    temp[triangle_pairs][0] =temp_triangle_node_C;
                    temp[triangle_pairs][1] =temp_triangle_node_A;
                    temp[triangle_pairs][2] =temp_triangle_node_B;
                    temp[triangle_pairs][3] =neighbor;
                }
                triangle_pairs++;
            }
            neighbor_indicator=0;
        }
    }
    
    for (int i=0; i<Membrane_triangle_list.size(); i++) {
        if (membrane_triangle_pair_list[i].size()!=3) {
            cout<<"\nThere is an error in the 'Membrane_Triangle_Pair_Identifier' function. This error indicates that there is a membrane triangle that has more/less than 3 triangle neighbors."<<endl;
//            cout<<"membrane_triangle_pair_list["<<i<<"].size()="<<membrane_triangle_pair_list[i].size()<<endl;
            cout<<"triangle "<<i<<" neighbors:\n";
            for (int k=0; k<membrane_triangle_pair_list[i].size(); k++) {
                cout<<membrane_triangle_pair_list[i][k]<<"\n";
//                for (int k2=0; k2<3; k2++) {
////                    Membrane_triangle_list[membrane_triangle_pair_list[i][k]][k2]<<"\t";
//                }
            }
        }
    }
    
    for(int i=0;i<2*Membrane_num_of_Triangle_Pairs;i++)//saving temp in temp2
    {
        for(int j=0;j<4;j++)
        {
            temp2[i][j]=temp[i][j];
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)// sorting temp
    {
        if( temp[abc][0]  > temp[abc][1]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][1] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][1] =temp_triangle_node_A;
        }
        if( temp[abc][0]  > temp[abc][2]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][2] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][2] =temp_triangle_node_A;
        }
        if( temp[abc][0]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
        
        if( temp[abc][1]  > temp[abc][2]    )
        {
            temp_triangle_node_A=temp[abc][1] ;
            temp_triangle_node_B=temp[abc][2] ;
            temp[abc][1] =temp_triangle_node_B;
            temp[abc][2] =temp_triangle_node_A;
        }
        
        if( temp[abc][1]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][1] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][1] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
        
        if( temp[abc][2]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][2] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][2] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if(temp[abc][0] !=-1)
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][1] ;
            temp_triangle_node_C=temp[abc][2] ;
            temp_triangle_node_D=temp[abc][3] ;
            for(int cab=0;cab<2*Membrane_num_of_Triangle_Pairs;cab++)
            {
                if( temp_triangle_node_A==temp[cab][0]  &   temp_triangle_node_B==temp[cab][1]    &   temp_triangle_node_C==temp[cab][2]    &   temp_triangle_node_D==temp[cab][3]  & abc!=cab  )
                {
                    temp[cab][0]  =-1;
                }
                
            }
            
        }
        
    }
    
    
    int temp_int=0;
	
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if( temp[abc][0] !=-1)
        {
            Membrane_Triangle_Pair_Nodes[temp_int][0]=temp2[abc][0] ;
            Membrane_Triangle_Pair_Nodes[temp_int][1]=temp2[abc][1] ;
            Membrane_Triangle_Pair_Nodes[temp_int][2]=temp2[abc][2] ;
            Membrane_Triangle_Pair_Nodes[temp_int][3]=temp2[abc][3] ;
			
			Membrane_Edges[temp_int][0]=temp2[abc][0];
			Membrane_Edges[temp_int][1]=temp2[abc][1];
			
            
            temp_int++;
        }
    }
}
