#include <iostream>
#include <cmath>
#include "cblas.h"
using namespace std;

// void zeros(int, float*);

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
    const int nelem_y = 5;    // Number of elements in y-direction
    const int nelem_x = 10;   // Number of elements in the x-direction

    int nnode_elem = 4 ;                // number of nodes in each element
    //float nNodeDof[4] = {1, 1, 1, 1} ;    // number of DoF per node (1 = Temperature only)
    // total number of DoF given by the sum of the nodal DoF
    //int neDof = cblas_sasum(4, nNodeDof, 1);    // CBLAS routine to sum elements of array nNodeDof

    //----- Element and node numbers -----------------
    int nelem = nelem_x * nelem_y;              // Total number of elements
    int nnode = (nelem_x + 1) * (nelem_y + 1);   // Total number of nodes




    // ------------------------
    // ** NODAL COORDINATES **
    // ------------------------

    //----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  // Beam height as a function of x

    // initialise Y with zeros using zeros function - not required, however kept for consideration
    float Y[(nelem_x+1) * (nelem_y+1)]; // ** replace Y index with nNode
    // float* Y_pointer = &(Y[0]);
    // zeros((nelem_x+1) * (nelem_y+1), Y_pointer);

    // create top and bottom position for each x and fill Y with h positions.
    // Y is stored as a column major matrix and hence is filled column by column.

    // initialise vector x to store nodal coordinates
    float x[nelem_x];
    // linear spacing of node coordinates
    for (int i = 0; i <= nelem_x; i++){

        x[i] = L / nelem_x * i;
    }

    // initialise h as the plate upper edge boundary (mirrored for lower boundary)
    float h[nelem_x];
    // calculate h for each point in x
    for (int i = 0; i <= nelem_x; i++) {

        h[i] = a * x[i] * x[i] + b * x[i] + h1;

    }




    for (int colnr = 0; colnr < nelem_x +1; colnr++){
        for (int rownr = 0; rownr < nelem_y + 1; rownr++){

            Y[colnr * nelem_x + rownr] = -h[colnr] / 2 + h[colnr] / nelem_y * rownr ;

        }
    }

    // create Coordinate matrix with two rows (x,y) for each node.
    float Coord[nnode * 2];

    // loop through each x and y value, producing coordinates. Stored as column major {x(:), y(:)}
    for (int i = 0; i < nelem_x+1; i++){
        for(int j = 0; j < nelem_y+1; j++){

            Coord[j + i*nelem_y] = x[i];                // first column : x
            Coord[(j + i*nelem_y) + nnode] = Y[j];      // second column : y
        }
    }

    // ----- Calculation of topology matrix NodeTopo -----------

    int NodeTopo[(nelem_x+1) * (nelem_y+1)];
    // since matrices are stored as column major, assignment of node numbers is simple
    for (int i = 0; i < nelem_x+1; i++){
        for (int j = 0; j < nelem_y+1; j++){

            NodeTopo[j + (nelem_y+1)*i] = j + (nelem_y+1)*i;

        }
    }

    // ----- Calculation of topology matrix ElemNode -----------

    // initialise array for topology matrix- 5 columns: 1x element number, 4x corners for a quadrilateral element
    int ElemNode[nelem*5];
    int elemnr = 0;

    for (int colnr = 0; colnr < nelem_x; colnr++){
        for (int rownr = 0; rownr < nelem_y; rownr++){

            ElemNode[elemnr] = elemnr;                                                                      // column 0: element number
            ElemNode[elemnr + (nelem_x * nelem_y) ] = NodeTopo[rownr + colnr * (nelem_y+1)];                 // column 1: top left - NodeTopo(row, column)
            ElemNode[elemnr + (nelem_x * nelem_y) * 2] = NodeTopo[rownr + colnr * (nelem_y+1) + nelem_y+1];             // column 2: top right - NodeTopo(row+1, column)
            ElemNode[elemnr + (nelem_x * nelem_y) * 3] = NodeTopo[rownr + colnr * (nelem_y+1) + nelem_y+1 + 1];           // column 3: bottom right - NodeTopo(row+1, column+1)
            ElemNode[elemnr + (nelem_x * nelem_y) * 4] = NodeTopo[rownr + colnr * (nelem_y+1) + 1];           // column 4: bottom left - NodeTopo(row+1, column)
            elemnr++;

        }

    }
    






    return 0;
}
/*
void zeros(int n, float* A){


    for (int i = 0; i < n ; i++){

        A[i] = 3;
        cout << A[i];

    }
}
 */
