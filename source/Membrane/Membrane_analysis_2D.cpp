#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>

using namespace std;

void Membrane::calculate_freqs(ArgStruct_Analysis args){
    
    int q_max = args.q_max;
    if(cn_2D_avg.size() != q_max+1){
//        Un2D_avg.clear();
        cn_2D_avg.clear();
        cn2_2D_avg.clear();
//        Un2D_avg.resize(q_max+1, 0);
        cn_2D_avg.resize(q_max+1, 0);
        cn2_2D_avg.resize(q_max+1, 0);
        if (!generalParameters.Testmode) {
            cout<<"cleared uq\n";
        }
    }
    
    
    if(args.analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    if (ringNodeList.size()==0) {
        get_ring(args);
    }
    
    
    
    
    
    calculate_contourSegmentLength();
    calculate_contourRadius();
    vector< double > cn_frame = calculate_2D_amplitudes(args.q_max);
    
    for (int i=0; i<=q_max; i++) {
        cn_2D_avg[i]+=cn_frame[i];
        cn2_2D_avg[i]+=cn_frame[i]*cn_frame[i];
    }
    
}



void Membrane::get_ring(ArgStruct_Analysis args){
    int n = int (sqrt(Num_of_Nodes))*2;
    double dtheta = M_PI/n;
    double theta0 = M_PI/2;
    double thresh_min = theta0 - dtheta;
    double thresh_max = theta0 + dtheta;
    
    vector<int> chainlist(Num_of_Nodes,0);
    vector<vector<double> > indPhi;
    int count=0;
    for (int i=0; i<Num_of_Nodes; i++) {
        double theta = spherical_positions[i][1];
        if (theta < thresh_max && theta > thresh_min) {
            vector <double> temp;
            temp.push_back(i);
            temp.push_back(spherical_positions[i][2]+M_PI);
            indPhi.push_back(temp);
            chainlist[i]=1;
            count++;
        }
        
    }
    //Sort node indecies based on their phi (from 0 to 2pi)
    sort(indPhi.begin(), indPhi.end(),[](const vector<double>& a, const vector<double>& b){
        return a[1] < b[1];
    });
    for (auto &i:indPhi) {
        ringNodeList.push_back(i[0]);
    }
    
    bool getPDBOutputForDoubleCheck=true;
    if (getPDBOutputForDoubleCheck) {
        cout<<"count "<<count<<endl;
        WriteMemPDBFrame(args,chainlist);
    }
    
    
}

void Membrane::WriteMemPDBFrame(ArgStruct_Analysis args,
                                vector<int> chainlist)
{
    string traj_name= args.output_filename[0];
    traj_name.erase(traj_name.end()-4, traj_name.end());
    traj_name+= "mem0.pdb";
    
    FILE* pFile;
    pFile = fopen (traj_name.c_str(),"a");
    fprintf(pFile,"MODEL     %d\n", 0);
    fprintf(pFile,"REMARK 250 time=%.3f ps;\n",
            0);
    char chain[]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
    for (int n=0; n<Num_of_Nodes; ++n){
        
        fprintf(pFile,"ATOM  %5d %4s ETH %c%4.0f    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                "mem0",
                chain[chainlist[n]],
                double(n+1),
                Node_Position[n][0],
                Node_Position[n][1],
                Node_Position[n][2],
                0,
                0);//,
        //                atoms[n].symbol);
    }


    fprintf(pFile,"ENDMDL\n");
    fclose (pFile);
}

void Membrane::calculate_contourSegmentLength(void){
    contourSegmentLength.clear();
    
    
    for (int i=0; i<ringNodeList.size()-1; i++) {
        double xi   = Node_Position[ringNodeList[i]][0];
        double yi   = Node_Position[ringNodeList[i]][1];
        double xip1 = Node_Position[ringNodeList[i+1]][0];
        double yip1 = Node_Position[ringNodeList[i+1]][1];
        contourSegmentLength.push_back( sqrt( (xip1-xi)*(xip1-xi) + (yip1-yi)*(yip1-yi) ) );
    }
    int N = int(ringNodeList.size()-1);
    double x0   = Node_Position[ringNodeList[0]][0];
    double y0   = Node_Position[ringNodeList[0]][1];
    double xN = Node_Position[ringNodeList[N]][0];
    double yN = Node_Position[ringNodeList[N]][1];
    contourSegmentLength.push_back( sqrt( (xN-x0)*(xN-x0) + (yN-y0)*(yN-y0) ) );
}

void Membrane::calculate_contourRadius(void){
    contourRadius=0;
    for (int i=0; i<ringNodeList.size()-1; i++) {
        double ri     = spherical_positions[ringNodeList[i]][0];
        double rip1   = spherical_positions[ringNodeList[i+1]][0];
        double phii   = spherical_positions[ringNodeList[i]][2];
        double phiip1 = spherical_positions[ringNodeList[i+1]][2];
        contourRadius += ( ri+rip1 )*( phiip1-phii );
        cout<<( phiip1-phii )<<endl;
    }
    int N = int(ringNodeList.size()-1);
    double r0     = spherical_positions[ringNodeList[0]][0];
    double rN   = spherical_positions[ringNodeList[N]][0];
    double phi0   = spherical_positions[ringNodeList[0]][2];
    double phiN = spherical_positions[ringNodeList[N]][2];
    contourRadius += ( r0+rN )*( 2*M_PI-phiN+phi0 );
    cout<<( 2*M_PI-phiN+phi0 )<<endl;
    contourRadius *= 1/(4*M_PI);
    cout<<"contourRadius "<<contourRadius<<endl;
}

vector<double> Membrane::calculate_2D_amplitudes(int q_max){
    vector<double> a_n(q_max+1,0);
    vector<double> b_n(q_max+1,0);
    vector<double> c_n(q_max+1,0);
    
    
    for (int q=0; q<=q_max; q++) {
        for (int i=0; i<ringNodeList.size()-1; i++) {
            double ri     = spherical_positions[ringNodeList[i]][0];
            double rip1   = spherical_positions[ringNodeList[i+1]][0];
            double phii   = spherical_positions[ringNodeList[i]][2];
            double phiip1 = spherical_positions[ringNodeList[i+1]][2];
            
            a_n[q]+= (phiip1-phii)*( ri*cos(q*phii) + rip1*cos(q*phiip1) );
            b_n[q]+= (phiip1-phii)*( ri*sin(q*phii) + rip1*sin(q*phiip1) );
        }
        
        int N = int(ringNodeList.size()-1);
        double r0     = spherical_positions[ringNodeList[0]][0];
        double rN   = spherical_positions[ringNodeList[N]][0];
        double phi0   = spherical_positions[ringNodeList[0]][2];
        double phiN = spherical_positions[ringNodeList[N]][2];
        
        a_n[q]+= (phiN-phi0)*( r0*cos(q*phi0) + rN*cos(q*phiN) );
        b_n[q]+= (phiN-phi0)*( r0*sin(q*phi0) + rN*sin(q*phiN) );
        
        a_n[q]/= 2*M_PI*contourRadius;
        b_n[q]/= 2*M_PI*contourRadius;
        
        c_n[q] = sqrt( a_n[q]*a_n[q] + b_n[q]*b_n[q]);
        
    }
    
    return c_n;
}

void Membrane::write_un_uq(ArgStruct_Analysis args, int file_index){
    
    double numOfFrames = args.framelimits_beg-args.framelimits_end;
    vector<double> Un2D_avg;
    vector<double> Uq2D_avg;
    Uq2D_avg.resize(args.q_max+1,0);
    Un2D_avg.resize(args.q_max+1,0);
    for (int i=0; i<=args.q_max; i++) {
        cn_2D_avg[i]/=numOfFrames;
        cn2_2D_avg[i]/=numOfFrames;
        Un2D_avg[i] = contourRadius*contourRadius*contourRadius*0.5*M_PI*(cn2_2D_avg[i]-cn_2D_avg[i]*cn_2D_avg[i]);
        Uq2D_avg[i] = (cn2_2D_avg[i]-cn_2D_avg[i]*cn_2D_avg[i]);
    }
    
    std::ofstream wdata;
    wdata.open(args.output_filename[file_index].c_str(), std::ios::app);
    wdata<<"#Average contour radius\n";
    wdata<<"#<|Un|^2> n=q*R\n";
    wdata<<"#<|Uq|^2>\n";
    
    wdata<<contourRadius<<endl;
    for (int n=0; n<args.q_max+1; n++) {
        wdata<<Un2D_avg[n]<<"\t";
    }
    wdata<<"\n";
    for (int q=0; q<args.q_max+1; q++) {
        wdata<<Uq2D_avg[q]<<"\t";
    }
    wdata<<endl;
}
