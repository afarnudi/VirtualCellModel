#include "Membrane.h"
#include "General_functions.hpp"
#include <cstdlib>



#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>


using namespace std::complex_literals;
using namespace std;

void Membrane::rotate_coordinates(double theta, double phi){
    vector<vector<double> > Rotationy;
    vector<vector<double> > Rotationz;
    
    Rotationy.resize(3);
    Rotationz.resize(3);
    for (int i=0; i<3; i++) {
        Rotationy[i].resize(3);
        Rotationz[i].resize(3);
    }
    if (theta != 0) {
        Rotationy[0][0] = cos(theta);
        Rotationy[0][1] = 0;
        Rotationy[0][2] = sin(theta);
        
        Rotationy[1][0] = 0;
        Rotationy[1][1] = 1;
        Rotationy[1][2] = 0;
        
        Rotationy[2][0] = -sin(theta);
        Rotationy[2][1] = 0;
        Rotationy[2][2] = cos(theta);
    }
    
    if (phi !=0) {
        Rotationz[0][0] = cos(phi);
        Rotationz[0][1] = -sin(phi);
        Rotationz[0][2] = 0;
        
        Rotationz[1][0] = sin(phi);
        Rotationz[1][1] = cos(phi);
        Rotationz[1][2] = 0;
        
        Rotationz[2][0] = 0;
        Rotationz[2][1] = 0;
        Rotationz[2][2] = 1;
    }
    
    
    for (int i=0; i<Num_of_Nodes; i++) {
        double x = Node_Position[i][0];
        double y = Node_Position[i][1];
        double z = Node_Position[i][2];
        
        if (phi != 0) {
            Node_Position[i][0] = x*Rotationz[0][0] + y*Rotationz[0][1] + z*Rotationz[0][2];
            Node_Position[i][1] = x*Rotationz[1][0] + y*Rotationz[1][1] + z*Rotationz[1][2];
            Node_Position[i][2] = x*Rotationz[2][0] + y*Rotationz[2][1] + z*Rotationz[2][2];
        }
        
        if (theta !=0) {
            x = Node_Position[i][0];
            y = Node_Position[i][1];
            z = Node_Position[i][2];
            
            Node_Position[i][0] = x*Rotationy[0][0] + y*Rotationy[0][1] + z*Rotationy[0][2];
            Node_Position[i][1] = x*Rotationy[1][0] + y*Rotationy[1][1] + z*Rotationy[1][2];
            Node_Position[i][2] = x*Rotationy[2][0] + y*Rotationy[2][1] + z*Rotationy[2][2];
        }
        
        
        
        
    }
    
    
    
}


void Membrane::update_spherical_positions(){
    spherical_positions.clear();
    spherical_positions.resize(Num_of_Nodes);
    
    for (int i=0; i<Num_of_Nodes; i++) {
        spherical_positions[i].resize(3,0);
        spherical_positions[i] = convert_cartesian_to_spherical(Node_Position[i][0],                                                                       Node_Position[i][1],
                                                                Node_Position[i][2]);
    }
}


void Membrane::rotate_particle_to_axes(void){
    
    //find the angles of the node that used to be on the z axis
    vector<double> r_theta_phi;
    r_theta_phi = convert_cartesian_to_spherical(Node_Position[z_node_index][0],
                                                 Node_Position[z_node_index][1],
                                                 Node_Position[z_node_index][2]);
    
    double theta0 = r_theta_phi[1];
    double phi0   = r_theta_phi[2];
    
    rotate_coordinates(-theta0, -phi0);
    //rotate the particles so the node that was on the z axis at the beginning of the run will lie on the z axis
    r_theta_phi = convert_cartesian_to_spherical(Node_Position[y_node_index][0],
                                                 Node_Position[y_node_index][1],
                                                 Node_Position[y_node_index][2]);
    
    phi0   = r_theta_phi[2];
    //rotate the particles so the node that was on the y axis at the beginning of the run will lie on the y axis
    rotate_coordinates(0, M_PI-phi0);
}


