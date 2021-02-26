#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>
#include <algorithm>

using namespace std;

void Membrane::calculate_freqs_projection(ArgStruct_Analysis args){
    int q_max = args.q_max;
    
    if(args.analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    int numpoints = 2*sqrt(Num_of_Nodes);
    deltaphi.clear();
    double phiinc=2*M_PI/numpoints;
    for (int i=0; i<numpoints; i++) {
        deltaphi.push_back(i*phiinc);
    }
    get_ring_from_projection(args, numpoints);
    contourRadius=0;
    for (int i=0; i<numpoints; i++) {
        contourRadius+=contourradii[i];
    }
    contourRadius/=numpoints;
    
    calculateH_2D_amplitudes_projection(args.q_max, numpoints);
    calculate_2D_amplitudes_Alexandra_projection(args.q_max, numpoints);
}

void Membrane::calculate_freqs(ArgStruct_Analysis args){
    
    int q_max = args.q_max;
   
    
    if(args.analysis_averaging_option == 1){
        double phi   = ((double) rand() / (RAND_MAX))*2*M_PI;
        double theta = ((double) rand() / (RAND_MAX))*M_PI;
        
        rotate_coordinates(theta, phi);
        update_spherical_positions();
        
    }
    if (ringNodeList.size()==0) {
        get_ring(args);
    }
    
    calculate_contourRadius();
    calculateH_2D_amplitudes(args.q_max);
    calculate_2D_amplitudes_Alexandra(args.q_max);
    
    bool getPDBOutputForDoubleCheck=false;
    if (getPDBOutputForDoubleCheck) {
        //        cout<<"count "<<count<<endl;
        WriteMemPDBFrame(args,chainlist);
    }
    
}


#include <iomanip>


void Membrane::calculateH_2D_amplitudes(int q_max){
    calculate_deltaphi_Alexandra();
    
    vector<double> a_q(q_max+1,0);
    vector<double> b_q(q_max+1,0);
    vector<double> c_q(q_max+1,0);
    vector<double> c_q2(q_max+1,0);
    
    for (int q=0; q<=q_max; q++) {
        
        for (int i=0; i<ringNodeList.size(); i++) {
            double ri     = (spherical_positions[ringNodeList[i]][0]-contourRadius);//contourRadius;
            double phii   = spherical_positions[ringNodeList[i]][2];
            
            //Note that for the integral discretisation the dphi at point i is 0.5*deltaphi, delta phi = phi_i+1 - phi_i-1
            a_q[q]+= ri*cos(q*phii)*0.5*deltaphi[i];
            b_q[q]+= ri*sin(q*phii)*0.5*deltaphi[i];
        }

        a_q[q]*= 1./(2*M_PI);
        b_q[q]*= 1./(2*M_PI);
        
    }
    Aq_H.push_back(a_q);
    Bq_H.push_back(b_q);
    
}

void Membrane::calculateH_2D_amplitudes_projection(int q_max, int numpoints){

    vector<double> a_q(q_max+1,0);
    vector<double> b_q(q_max+1,0);
    vector<double> c_q(q_max+1,0);
    vector<double> c_q2(q_max+1,0);
    
    double deltaphi_projection = 2*M_PI/numpoints;
    
    for (int q=0; q<=q_max; q++) {
        
        for (int i=0; i<contourradii.size(); i++) {
            double ri     = (contourradii[i]-contourRadius);//contourRadius;
            double phii   = i*deltaphi_projection;
            
            //Note that for the integral discretisation the dphi at point i is 0.5*deltaphi, delta phi = phi_i+1 - phi_i-1
            a_q[q]+= ri*cos(q*phii);
            b_q[q]+= ri*sin(q*phii);
        }
        a_q[q]*= deltaphi_projection/(2*M_PI);
        b_q[q]*= deltaphi_projection/(2*M_PI);

    }
    Aq_H_proj.push_back(a_q);
    Bq_H_proj.push_back(b_q);
    
}

void Membrane::calculate_2D_amplitudes_Alexandra(int q_max){
    calculate_deltaphi_Alexandra();
    vector<double> a_q(q_max+1,0);
    vector<double> b_q(q_max+1,0);
    
    
    for (int q=0; q<=q_max; q++) {
        
        for (int i=0; i<ringNodeList.size(); i++) {
            double ri     = (spherical_positions[ringNodeList[i]][0]-contourRadius)/contourRadius;
            double phii   = spherical_positions[ringNodeList[i]][2];
            
            //Note that for the integral discretisation the dphi at point i is 0.5*deltaphi, delta phi = phi_i+1 - phi_i-1
            a_q[q]+= ri*cos(q*phii)*0.5*deltaphi[i];
            b_q[q]+= ri*sin(q*phii)*0.5*deltaphi[i];
        }
        a_q[q]/= (2*M_PI);
        b_q[q]/= (2*M_PI);
        
    }
    aq_alexandra.push_back(a_q);
    bq_alexandra.push_back(b_q);
    
}

void Membrane::calculate_2D_amplitudes_Alexandra_projection(int q_max, int numpoints){
    
    vector<double> a_q(q_max+1,0);
    vector<double> b_q(q_max+1,0);
    
    double deltaphi_projection = 2*M_PI/numpoints;
    
    for (int q=0; q<=q_max; q++) {
        
        for (int i=0; i<contourradii.size(); i++) {
            double ri     = (contourradii[i]-contourRadius)/contourRadius;
            double phii   = i*deltaphi_projection;
            
            //Note that for the integral discretisation the dphi at point i is 0.5*deltaphi, delta phi = phi_i+1 - phi_i-1
            a_q[q]+= ri*cos(q*phii);
            b_q[q]+= ri*sin(q*phii);
        }
        a_q[q]*= deltaphi_projection/(2*M_PI);
        b_q[q]*= deltaphi_projection/(2*M_PI);
        
    }
    aq_alexandra_proj.push_back(a_q);
    bq_alexandra_proj.push_back(b_q);
    
}

void Membrane::updatepos(double scale){
    for (int i=0; i<Num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            Node_Position[i][j]*=scale;
        }
    }
}

void Membrane::get_ring_from_projection(ArgStruct_Analysis args, int numpoints){
    if(args.analysis_averaging_option == 2){
        if(args.z_node == -1 || args.zy_node == -1){
            cout<<"No nodes specified for aligning.\n the nodes must be specfied with the '--align_axes flag z,y' where z is the node index to be aligned with the  z axis (when the com is in the origin (0,0,0) )and y is the index of the node that will then be rotated to lie on the zy plane.\n";
            exit(EXIT_FAILURE);
        }
        rotate_particle_to_axes(args);
    }
    update_COM_position();
    set_com_to_zero();
    update_spherical_positions();
    
    
    contourradii.clear();
    contourradii.resize(numpoints,0);
    
    double phiinc = 2*M_PI/numpoints;
    
    for (int i=0; i<Num_of_Nodes; i++) {
        double projection_radius = sqrt(Node_Position[i][0]*Node_Position[i][0]+Node_Position[i][1]*Node_Position[i][1]);
        int radius_index = (M_PI+spherical_positions[i][2])/phiinc;
        if (contourradii[radius_index]<projection_radius) {
            contourradii[radius_index]=projection_radius;
        }
        
    }
    for (int i=0; i<numpoints; i++) {
        
        if (contourradii[i]==0) {
            string errorMessage = TWARN;
            errorMessage+="Error: get_ring_from_projection: got a zero radius. Decrease the num of points you detect on the contour. Current num of points is set to "+to_string(numpoints)+"."; errorMessage+=TRESET;
            errorMessage+="\n";
            throw std::runtime_error(errorMessage);
        }
    }
    
    
}

void Membrane::get_ring(ArgStruct_Analysis args){
    if(args.analysis_averaging_option == 2){
        if(args.z_node == -1 || args.zy_node == -1){
            cout<<"No nodes specified for aligning.\n the nodes must be specfied with the '--align_axes flag z,y' where z is the node index to be aligned with the  z axis (when the com is in the origin (0,0,0) )and y is the index of the node that will then be rotated to lie on the zy plane.\n";
            exit(EXIT_FAILURE);
        }
        rotate_particle_to_axes(args);
    }
    update_COM_position();
    set_com_to_zero();
    update_spherical_positions();
    
    
    
    
    int n = int (sqrt(Num_of_Nodes))*2;
    double dtheta_persision=0.9;
    double dtheta = M_PI/n;
    double theta0 = M_PI/2;
    double thresh_min;// = theta0 - dtheta*dtheta_persision;
    double thresh_max;// = theta0 + dtheta*dtheta_persision;
    
    vector<vector<double> > indPhi;
    dtheta_persision+=0.1;
//    cout<<"here, persision "<<dtheta_persision<<" indphi: "<<indPhi.size()<<endl;
    thresh_min = theta0 - dtheta*dtheta_persision;
    thresh_max = theta0 + dtheta*dtheta_persision;
    chainlist.clear();
    chainlist.resize(Num_of_Nodes,0);
    indPhi.clear();
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
    if (indPhi.size()<0.78*n) {
        while (indPhi.size()<n) {
            
            dtheta_persision+=0.1;
            cout<<"here, persision "<<dtheta_persision<<" indphi: "<<indPhi.size()<<endl;
            thresh_min = theta0 - dtheta*dtheta_persision;
            thresh_max = theta0 + dtheta*dtheta_persision;
            chainlist.clear();
            chainlist.resize(Num_of_Nodes,0);
            indPhi.clear();
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
            if (dtheta_persision>5) {
                string errorMessage = TWARN;
                errorMessage+="Error: get_ring: cannot detect a ring of nodes on the equator."; errorMessage+=TRESET;
                errorMessage+="\n";
                throw std::runtime_error(errorMessage);
            }
        }
    
    }
    
    //Sort node indecies based on their phi (from 0 to 2pi)
    std::sort(indPhi.begin(), indPhi.end(),[](const vector<double>& a, const vector<double>& b){
        return a[1] < b[1];
    });
    for (auto &i:indPhi) {
        ringNodeList.push_back(i[0]);
//        cout<<i[0]<<" ";//<<i[1]<<endl;
    }
    cout<<endl;
    bool getPDBOutputForDoubleCheck=true;
    if (getPDBOutputForDoubleCheck) {
        //        cout<<"count "<<count<<endl;
        WriteMemPDBFrame(args,chainlist);
    }
    cout<<"\nIdentified "<<ringNodeList.size()<<" nodes on the membrane contour\n"<<endl;
//    ringNodeList.clear();
}

void Membrane::WriteMemPDBFrame(ArgStruct_Analysis args,
                                vector<int> list)
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
    int index=0;
    for (int n=0; n<Num_of_Nodes; ++n){
        
        fprintf(pFile,"ATOM  %5d %4s ETH %c   0    %8.3f%8.3f%8.3f%6.2f%6.1f\n",
                n+1,
                "mem0",
                chain[list[n]],
//                index,
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

void Membrane::calculate_deltaphi(void){
    deltaphi.clear();
    deltaphi.resize(ringNodeList.size(),0);
    
    for (int i=0; i<ringNodeList.size()-1; i++) {
        double phii   = spherical_positions[ringNodeList[i]][2];
        double phiip1 = spherical_positions[ringNodeList[i+1]][2];
        deltaphi[i] = ( phiip1-phii );
        if (deltaphi[i]<0) {

            if (deltaphi[i] < -M_PI) {
                deltaphi[i] += 2*M_PI;
            } else {
                deltaphi[i]*=-1;
            }

        }
    }
    int N = int(ringNodeList.size()-1);
    
    double phi0   = spherical_positions[ringNodeList[0]][2];
    double phiN = spherical_positions[ringNodeList[N]][2];
    deltaphi[N] = (phi0-phiN);
    if (deltaphi[N]<0) {

        if (deltaphi[N] < -M_PI) {
            deltaphi[N] += 2*M_PI;
        } else {
            deltaphi[N]*=-1;
        }
        

    }
    
    for (auto &dphi:deltaphi) {
        if (dphi<0 || dphi>M_PI) {
            string errorMessage = TWARN;
            errorMessage+="Error: calculate_deltaphi for 2D analysis: "; errorMessage+=TRESET;
            errorMessage+="A delta phi calculation resulted in an unexpected value, got "; errorMessage+=TFILE;
            errorMessage+=to_string(dphi); errorMessage+=TRESET;
            errorMessage+="\n";
            throw std::runtime_error(errorMessage);
        }
    }
    
}

void Membrane::calculate_deltaphi_Alexandra(void){
    deltaphi.clear();
    deltaphi.resize(ringNodeList.size(),0);
    
    int N = int(ringNodeList.size()-1);
    
    double phiN = spherical_positions[ringNodeList[N]][2];
    double phi1 = spherical_positions[ringNodeList[1]][2];
    deltaphi[0] = (phi1-phiN);
    
    if (deltaphi[0]<0) {
        if (deltaphi[0] < -M_PI) {
            deltaphi[0] += 2*M_PI;
        } else {
            deltaphi[0]*=-1;
        }
    }
    
    for (int i=1; i<ringNodeList.size()-1; i++) {
        double phiim1 = spherical_positions[ringNodeList[i-1]][2];
        double phiip1 = spherical_positions[ringNodeList[i+1]][2];
        deltaphi[i] = ( phiip1-phiim1 );
        
        if (deltaphi[i]<0) {
            if (deltaphi[i] < -M_PI) {
                deltaphi[i] += 2*M_PI;
            } else {
                deltaphi[i]*=-1;
            }
        }
    }
    
    
    double phiNm1 = spherical_positions[ringNodeList[N-1]][2];
    double phi0   = spherical_positions[ringNodeList[0]][2];
    deltaphi[N] = (phi0-phiNm1);
    if (deltaphi[N]<0) {
        if (deltaphi[N] < -M_PI) {
            deltaphi[N] += 2*M_PI;
        } else {
            deltaphi[N]*=-1;
        }
    }
    
    for (auto &dphi:deltaphi) {
        if (dphi<0 || dphi>M_PI) {
            string errorMessage = TWARN;
            errorMessage+="Error: calculate_deltaphi_alexandra for 2D analysis: "; errorMessage+=TRESET;
            errorMessage+="A delta phi calculation resulted in an unexpected value, got "; errorMessage+=TFILE;
            errorMessage+=dphi; errorMessage+=TRESET;
            errorMessage+="\n";
            throw std::runtime_error(errorMessage);
        }
    }
}


void Membrane::calculate_contourRadius(void){
    calculate_deltaphi();
    contourRadius=0;
    for (int i=0; i<ringNodeList.size()-1; i++) {
        double ri     = spherical_positions[ringNodeList[i]][0];
        double rip1   = spherical_positions[ringNodeList[i+1]][0];
        contourRadius += ( ri+rip1 )*deltaphi[i];
        
    }
    int N = int(ringNodeList.size()-1);
    double r0     = spherical_positions[ringNodeList[0]][0];
    double rN   = spherical_positions[ringNodeList[N]][0];
    contourRadius += ( r0+rN )*deltaphi[N];
    
    contourRadius /= 4*M_PI;
}



void Membrane::write_un_uq(ArgStruct_Analysis args, int file_index){
    
    double numOfFrames = args.framelimits_end-args.framelimits_beg;

    vector<double> Aq_avg;
    vector<double> Bq_avg;
    vector<double> Uq2_avg;
    Aq_avg.resize(args.q_max+1,0);
    Bq_avg.resize(args.q_max+1,0);
    Uq2_avg.resize(args.q_max+1,0);
    double M = numOfFrames;
    for (int i=0; i<M; i++) {
        for (int q=2; q<=args.q_max; q++) {
            Aq_avg[q] += Aq_H[i][q]/M;
            Bq_avg[q] += Bq_H[i][q]/M;
        }
    }
    
    for (int i=0; i<M; i++) {
        for (int q=0; q<=args.q_max; q++) {
            Uq2_avg[q] += ( ( Aq_H[i][q] - Aq_avg[q] )*( Aq_H[i][q] - Aq_avg[q] ) + ( Bq_H[i][q] - Bq_avg[q] )*( Bq_H[i][q] - Bq_avg[q] ) )/M;
            
        }
    }
    
    
    std::ofstream wdata;
    wdata.open(args.output_filename[file_index].c_str(), std::ios::app);
    wdata<<"#Average contour radius\n";
    wdata<<"#<|Un|^2> n=q*R\n";
    wdata<<"#<|Uq|^2> Standard Fourier transform \n";
    
    wdata<<contourRadius<<endl;
    for (int n=0; n<args.q_max+1; n++) {
//        double cq2_avg = 0.5*(Aq2_avg[n]+Bq2_avg[n]-Aq_avg[n]*Aq_avg[n]-Bq_avg[n]*Bq_avg[n]);
//        wdata<<cq2_avg<<"\t";
        wdata<<Uq2_avg[n]<<"\t";
    }
    
    
    
    M = args.framelimits_end-args.framelimits_beg;
    if( M != aq_alexandra.size()){
        cout<<"frame num = "<<M<<"  t = "<<aq_alexandra.size()<<endl;
        exit(0);
    }
    
    vector<double> aq_avg;
    vector<double> bq_avg;
    vector<double> uq2_avg;
    aq_avg.resize(args.q_max+1,0);
    bq_avg.resize(args.q_max+1,0);
    uq2_avg.resize(args.q_max+1,0);
    
    for (int i=0; i<M; i++) {
        for (int q=2; q<=args.q_max; q++) {
            aq_avg[q] += aq_alexandra[i][q]/M;
            bq_avg[q] += bq_alexandra[i][q]/M;
        }
    }
    
    for (int i=0; i<M; i++) {
        for (int q=0; q<=args.q_max; q++) {
            uq2_avg[q] += ( ( aq_alexandra[i][q] - aq_avg[q] )*( aq_alexandra[i][q] - aq_avg[q] ) + ( bq_alexandra[i][q] - bq_avg[q] )*( bq_alexandra[i][q] - bq_avg[q] ) )/M;
            
        }
    }
    
    wdata<<"\n";
    for (int q=0; q<args.q_max+1; q++) {
        wdata<<uq2_avg[q]<<"\t";
    }
    wdata<<endl;
}

void Membrane::write_un_uq_projection(ArgStruct_Analysis args, int file_index){
    
    double numOfFrames = args.framelimits_end-args.framelimits_beg;

    vector<double> Aq_avg;
    vector<double> Bq_avg;
    vector<double> Aq2_avg;
    vector<double> Bq2_avg;
    
    Aq_avg.resize(args.q_max+1,0);
    Bq_avg.resize(args.q_max+1,0);
    Aq2_avg.resize(args.q_max+1,0);
    Bq2_avg.resize(args.q_max+1,0);
    
    for (int t=0;t<numOfFrames; t++){
        for (int i=0; i<=args.q_max; i++) {
        
            Aq_avg[i]+=Aq_H_proj[t][i];
            Bq_avg[i]+=Bq_H_proj[t][i];
            Aq2_avg[i]+=Aq_H_proj[t][i]*Aq_H_proj[t][i];
            Bq2_avg[i]+=Bq_H_proj[t][i]*Bq_H_proj[t][i];
        }
        
    }
    
    for (int i=0; i<=args.q_max; i++) {
        Aq_avg[i]/=numOfFrames;
        Bq_avg[i]/=numOfFrames;
        Aq2_avg[i]/=numOfFrames;
        Bq2_avg[i]/=numOfFrames;
    }
    
    
    
    
    std::ofstream wdata;
    string projection_name = args.output_filename[file_index];
    
    projection_name.erase(projection_name.end()-7,projection_name.end());
    projection_name+="proj.txt";
    
    
    wdata.open(projection_name.c_str(), std::ios::app);
    wdata<<"#Average contour radius of the projected sphere onto the xyplane\n";
    wdata<<"#<|Un|^2> n=q*R\n";
    wdata<<"#<|Uq|^2> Standard Fourier transform \n";
    
    wdata<<contourRadius<<endl;
    
    for (int n=0; n<args.q_max+1; n++) {
        double cq2_avg = 0.5*(Aq2_avg[n]+Bq2_avg[n]-Aq_avg[n]*Aq_avg[n]-Bq_avg[n]*Bq_avg[n]);
        wdata<<cq2_avg<<"\t";
    }
    
    
    
    double M = args.framelimits_end-args.framelimits_beg;
    if( M != aq_alexandra_proj.size()){
        cout<<"frame num = "<<M<<"  t = "<<aq_alexandra_proj.size()<<endl;
        exit(0);
    }
    
    vector<double> aq_avg;
    vector<double> bq_avg;
    vector<double> uq2_avg;
    aq_avg.resize(args.q_max+1,0);
    bq_avg.resize(args.q_max+1,0);
    uq2_avg.resize(args.q_max+1,0);
    
    for (int i=0; i<M; i++) {
        for (int q=2; q<=args.q_max; q++) {
            aq_avg[q] += aq_alexandra_proj[i][q]/M;
            bq_avg[q] += bq_alexandra_proj[i][q]/M;
        }
    }
    
    for (int i=0; i<M; i++) {
        for (int q=0; q<=args.q_max; q++) {
            uq2_avg[q] += ( ( aq_alexandra_proj[i][q] - aq_avg[q] )*( aq_alexandra_proj[i][q] - aq_avg[q] ) + ( bq_alexandra_proj[i][q] - bq_avg[q] )*( bq_alexandra_proj[i][q] - bq_avg[q] ) )/M;
            
        }
    }
    
    wdata<<"\n";
    for (int q=0; q<args.q_max+1; q++) {
        wdata<<uq2_avg[q]<<"\t";
    }
    wdata<<endl;
}




#include <boost/math/special_functions/factorials.hpp>
void Membrane::write_uq_SH(ArgStruct_Analysis args, int file_index){
    
    double numOfFrames = args.framelimits_end-args.framelimits_beg;
    
    for (int ell=0; ell<args.ell_max+1; ell++) {
        for (int m=-ell; m<ell+1; m++) {
            
            ulm_avg[ell][m+ell]/=numOfFrames;
        }
    }
    vector<double> Uq2_avg;
    Uq2_avg.resize(args.q_max+1,0);
    double ell_max=args.ell_max;
    
    vector<vector<double> > ALP;
    ALP.resize(ell_max+1);
    for (int ell=0; ell<ell_max+1; ell++) {
        ALP[ell].resize(2*ell+1,0);
        for (int m=-ell; m<ell_max+1; m++) {
            ALP[ell][m+ell] = calc_assoc_legendre(ell,m,0);
            ALP[ell][m+ell]*= calc_assoc_legendre(ell,m,0);
            
            ALP[ell][m+ell]*= 2*ell+1;
            ALP[ell][m+ell]/= M_PI;
            
            ALP[ell][m+ell]*= boost::math::factorial<double>(ell-m);
            ALP[ell][m+ell]/= boost::math::factorial<double>(ell+m);
            
            
        }
    }
    
    for (int q=2; q<=args.q_max; q++) {
        for (int ell=q; ell<ell_max; q++) {
            Uq2_avg[q] += ALP[ell][q+ell]*ulm_avg[ell][q+ell];
        }
        
    }
    
    std::ofstream wdata;
    wdata.open(args.output_filename[file_index].c_str(), std::ios::app);
    wdata<<"#Average contour radius\n";
    wdata<<"#<|Uq|^2> calculated using the legendre polynomilas and ulms";
    wdata<<"#<|Uq|^2>\n";
    
    wdata<<0<<endl;
    for (int n=0; n<args.q_max+1; n++) {
        wdata<<Uq2_avg[n]<<"\t";
    }
    wdata<<"\n";
    for (int q=0; q<args.q_max+1; q++) {
        wdata<<Uq2_avg[q]<<"\t";
    }
    wdata<<endl;
    
}
