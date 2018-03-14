#include <iostream>
#include <cmath>
#include "cblas.h"
using namespace std;

int main() {

    // -----------------------
    // ** INITIALISE VALUES **
    // -----------------------

    // ----- Defining Materials -----------------
    float kx = 250.0;  // Thermal conductivity [W/mK]
    float ky = 250.0;  // Thermal conductivity [W/mK]
    float kxy = 0.0 ;  // Thermal conductivity [W/mK]


    // ----- Defining Section -----------------
    float th = .2;       // Thickness [m]

    // ----- Integration scheme -----------------
    int gaussorder = 2;


    // ----- Defining Geometry -----------------
    // h = a * x^2 + b * x + h1:    Beam height as a function of x
    float a = 0.0;          // constant in the polynomial describing the beam height
    float h1 = 1.0;         // [m] Height at beam left edge
    float h2 = 1.0 * h1;    // [m] Height at beam right edge
    float L = 2.0 * h1;     // [m] Beam length

    // calculate b constant in the polynomial describing the beam height
    float b = -a * L + (h2 - h1) / L;

    //----- Meshing Geometry -----------------
    int nelem_y = 8;    // Number of elements in y-direction
    int nelem_x = 15;   // Number of elements in the x-direction

    int nnode_elem = 4 ;                // number of nodes in each element
    double nNodeDof[4] = {1, 1, 1, 1} ;    // number of DoF per node (1 = Temperature only)
    // total number of DoF per element given by the sum the nodal DoF
    int neDof = cblas_dasum(4, nNodeDof, 1);    // dot product with an array of ones to sum elements of array nNodeDof
    cout << neDof;
    //----- Element and node numbers -----------------
    int nelem = nelem_x * nelem_y;              // Total number of elements
    int nnode = (nelem_x + 1) * (nelem_y + 1);   // Total number of nodes




    // ------------------------
    // ** NODAL COORDINATES **
    // ------------------------

    //----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  // Beam height as a function of x
//                                       x = np.linspace(0, L, nelem_x + 1)

//    h = a * x**2 + b * x + h1
//    Y = np.zeros((nelem_y + 1, nelem_x + 1))
//
//    for colnr in range(nelem_x + 1):
//    Y[:, colnr] = np.linspace(-h[colnr] / 2, h[colnr] / 2, nelem_y + 1) // create top and bottom position for each x and fill Y with h positions.
//
//            Coord = np.zeros((nnode, 2))  // Initialise coordinate of the nodes [x, y]



    return 0;
}