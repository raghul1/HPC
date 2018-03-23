#include <iostream>
#include <cmath>
#include "cblas.h"
#include "callfunctions.h"
#include <mpi.h>
#include <fstream>
#include <iomanip>
using namespace std;
#define F77NAME(x) x##_

extern "C" {
// LAPACK routine for solving systems of linear equations
void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B,
                    const int& ldb, int& info);
}

int main(int argc, char* argv[]) {

    /// **Takes command line input for a specified geometry, loading and boundary condition
    /// in order to solve the heat equation**.
    /** Receives command line input
     * in order to solve conductance equation. Using a specified number of quadrilateral elements, Gaussian integration is
     * used. The stiffness matrix is formed and later split into different sections. Finally, the output is written into a vtk file.*/

    ///@param case_no   case number (1, 2, 3) corresponding to the conditions given in the handout
    ///@param a         x^2 coefficient in polynomial describing height wrt x
    ///@param h1        height at beam left edge [m]
    ///@param h2        height at beam right edge [m]
    ///@param L         beam length [m]
    ///@param th        beam thickness [m]
    ///@param nelem_x   number of elements in x-direction
    ///@param nelem_y   number of elements in the y-direction
    ///@param T0        boundary condition temperature at boundary condition edge
    ///@param q         flux into edge
    ///@param kx, ky, kxy   thermal conductivity matrix values


    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank); // zero-based
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        cout << "Invalid communicator" << endl;
        return 1; }

    cout << "I am process " << rank + 1 << " of " << size << endl;

    // case_no, a, h1, h2, L, th, nelem_x, nelem_y, T0, q, kx, ky, kxy

    int case_no = atoi(argv[1]);
    // ----- Defining Geometry -----------------
    // h = a * x^2 + b * x + h1:    Beam height as a function of x
    float a = atof(argv[2]);          // constant in the polynomial describing the beam height
    float h1 = atof(argv[3]);         // [m] Height at beam left edge
    float h2 = atof(argv[4]);         // [m] Height at beam right edge
    float L = atof(argv[5]);         // [m] Beam length

    // ----- Defining section -----------------
    float th = atof(argv[6]);       // Thickness [m]

    // calculate b constant in the polynomial describing the beam height
    float b = -a * L + (h2 - h1) / L;
    const int nelem_y = atoi(argv[7]);    // Number of elements in y-direction
    const int nelem_x = atoi(argv[8]);   // Number of elements in the x-direction

    // Temperature and flux
    int T0 = atoi(argv[9]);
    int q = atoi(argv[10]);



    // ----- Defining material properties -----------------
    float kx = atof(argv[11]);  // Thermal conductivity [W/mK]
    float ky = atof(argv[12]);  // Thermal conductivity [W/mK]
    float kxy = atof(argv[13]); // Thermal conductivity [W/mK]


    // ----- Integration scheme -----------------
    int gaussorder = 2;

    int nnode_elem = 4;                // number of nodes in each element
    int neDof = 0;
    int nNodeDof[4] = {1, 1, 1, 1};    // number of DoF per node (1 = Temperature only)

    // total number of DoF given by the sum of the nodal DoF
    // int neDof = cblas_sasum(4, nNodeDof, 1);    // CBLAS routine to sum elements of array nNodeDof
    for (int i = 0; i < 4; i++) {
        neDof += nNodeDof[i];
    }

    //----- Element and node numbers -----------------
    int nelem = nelem_x * nelem_y;              // Total number of elements
    const int nnode = (nelem_x + 1) * (nelem_y + 1);   // Total number of nodes
    int nnode_x = nelem_x + 1;
    int nnode_y = nelem_y + 1;

    double D[4] = {kx, kxy, kxy, ky};                        // conductivity matrix D stored as column major

    // ** CHECKS **
    // check lengths given > 0m
    if (h1 <=0 || h2 <= 0 || L <=0 || th <= 0){
        cout << "Non-positive length entered. Please check dimensions." << endl;
        return 1;
    }
    // check number of elements in x & y > 0
    if (nelem_x <=0 || nelem_y <= 0){
        cout << "Non-positive number of elements entered. Please check element numbers." << endl;
        return 1;
    }

    // 2x2 matrix is positive if det > 0
    float detD = D[0] * D[3] - D[1] * D[2];
    if (detD <= 0){
        cout << "Matrix is not positive definite. Please check conductivity matrix values." << endl;
        return 1;
    }


    // -----------------------
    // **       CASES       **
    // -----------------------
    // initialisation for different cases

    int nFluxNodes;
    int nTempNodes;

    // set temperature and flux node numbers based on case number.
    // Case 1 & 2: tempnodes and fluxnodes are the same since conditions applied on opposite sides, but not case 3
    switch (case_no) {

            /// case 1: T_LEFT & Q_RIGHT
        case 1: {

            nFluxNodes = nnode_y;
            nTempNodes = nnode_y;
            break;

        };
            /// case 2: T_BOTTOM & Q_LEFT
        case 2: {

            nFluxNodes = nnode_x;
            nTempNodes = nnode_x;
            break;

        }

            /// case 3: T_LEFT & Q_BOTTOM
        case 3: {

            nFluxNodes = nnode_x;
            nTempNodes = nnode_y;
            break;

        }

    }

    // ------------------------
    // ** NODAL COORDINATES **
    // ------------------------

    //----- Calculation of Nodal coordinate matrix -> Coord ---------
    // h(x) = a * x**2 + b * x + h1  // Beam height as a function of x

    // initialise Y with zeros using zeros function - not required, however kept for consideration
    // float* Y_pointer = &(Y[0]);
    // zeros((nelem_x+1) * (nelem_y+1), Y_pointer);

    // create top and bottom position for each x and fill Y with h positions.
    // Y is stored as a column major matrix and hence is filled column by column.

    // initialise vector x to store nodal coordinates
    float x[nelem_x];
    float h[nelem_x];
    for (int i = 0; i < nnode_x; i++) {
        x[i] = L / nelem_x * i;                     // linear spacing of node coordinates
        h[i] = a * x[i] * x[i] + b * x[i] + h1;     // calculate h for each point in x
    }

    float Y[nnode_x * nnode_y];
    for (int colnr = 0; colnr < nnode_x; colnr++) {
        for (int rownr = 0; rownr < nnode_y; rownr++) {

            Y[colnr * (nelem_y + 1) + rownr] = (-h[colnr] / 2) + (h[colnr] / nelem_y) * rownr;
        }
    }

    // create coordinate matrix with two rows (x,y) for each node.
    // create topoology matrix defining position of numbered nodes
    double Coord[nnode * 2];
    int NodeTopo[nnode_x * nnode_y];
    for (int i = 0; i <
                    nnode_x; i++) {                  // loop through each x and y value, producing coordinates. Stored as column major x(:), y(:)
        for (int j = 0; j < nnode_y; j++) {

            Coord[j + i * nnode_y] = x[i];              // first column : x
            Coord[j + i * nnode_y + nnode] = Y[j];      // second column : y
            NodeTopo[j + i * nnode_y] = j + nnode_y * i;   // assignment of node numbers is trivial due to column major storage

        }
    }


    // ----- Calculation of topology matrix ElemNode -----------

    // initialise array for topology matrix- 5 columns: 1x element number, 4x corners for a quadrilateral element
    int ElemNode[nelem * 5];
    int elemnr = 0;

    for (int colnr = 0; colnr < nelem_x; colnr++) {
        for (int rownr = 0; rownr < nelem_y; rownr++) {

            ElemNode[elemnr] = elemnr;                                                        // column 0: element number
            ElemNode[elemnr + nelem] = NodeTopo[rownr + colnr *
                                                        nnode_y];                 // column 1: top left - NodeTopo(row, column)
            ElemNode[elemnr + nelem * 2] = NodeTopo[rownr + (colnr + 1) *
                                                            nnode_y];             // column 2: top right - NodeTopo(row+1, column)
            ElemNode[elemnr + nelem * 3] = NodeTopo[rownr + 1 + (colnr + 1) *
                                                                nnode_y];           // column 3: bottom right - NodeTopo(row+1, column+1)
            ElemNode[elemnr + nelem * 4] = NodeTopo[rownr + 1 + colnr *
                                                                nnode_y];               // column 4: bottom left - NodeTopo(row+1, column)
            elemnr++;

        }
    }



    // ----- Generate global DoF numbers -----------------
    // global DoF per node
    int globDof[nnode * 2];                     // globDof in form (global DoF, node number)
    zeros_int(nnode * 2, &(globDof[0]));        // zero matrix
    int temp_node = 0;                          // temporary node reassigned each loop

    for (int i = 0; i < nelem; i++) {
        for (int k = 0; k < 4; k++) {
            temp_node = ElemNode[k * nelem + i];
            // if the already existing ndof of the present node (global) is less than the present elements' ndof (local) then replace
            if (globDof[temp_node] < nNodeDof[k]) {
                globDof[temp_node] = nNodeDof[k];
            }
        }
    }

    // counting the global DoFs and inserting in globDof
    int nDof = 0;

    // for each node, add the number of Dof into nDof and insert node number into second row of globalDof
    for (int i = 0; i < nnode; i++) {
        for (int k = 0; k < globDof[i]; k++) {
            globDof[i + nnode] = nDof;
            nDof++;
        }
    }


    // ------------------------
    // **      K MATRIX     **
    // ------------------------

    // ----- Assembly of global stiffness matrix K ---------------

    double GP[2] = {-1 / sqrt(3.0), 1 / sqrt(3.0)};         // Gauss points
    int W[2] = {1, 1};                                      // Weights


    double K[nnode * nnode];                                 // global stiffness matrix K
    zeros_double(nnode * nnode, &(K[0]));                    // due to iterative addition, K must be zeroed to prevent unpredictable behaviour
    int eNodes[nnode_elem];                                 // element node numbers
    double eCoord[nnode_elem * 2];
    int gDof[nnode_elem];                                   // local element DoFs
    double Ke[nnode_elem * nnode_elem];                     // local element stiffness matrix
    double N[nnode_elem];                                   // shape function matrix
    double GN[nnode_elem * gaussorder];                     // derivative (gradient) of shape functions matrix
    double J[gaussorder * gaussorder];                      // Jacobian
    double detJ;                                            // Jacobian determinate
    double invJ[gaussorder * gaussorder];                   // Inverse Jacobian
    double B[nnode_elem * gaussorder];                      // Strain interpolation matrix
    double eta;
    double xi;
    double BtD[8];                                          // B' * D
    double alpha;                                           //  coefficient th * DetJ * W[i] * W[j]


    // Data for each node j in element i
    for (int i = 0; i < nelem; i++) {
        zeros_double(nnode_elem * nnode_elem,
                     &(Ke[0]));                  // remember to zero Ke after each element iteration

        for (int j = 0; j < nnode_elem; j++) {
            eNodes[j] = ElemNode[nelem * (j + 1) + i];                  // Node numbers
            eCoord[j] = Coord[eNodes[j]];                               // x (at j) and y (at j+nnode_elem)
            eCoord[nnode_elem + j] = Coord[eNodes[j] + nnode];
            gDof[j] = globDof[eNodes[j] + nNodeDof[j] * nnode];         // global Dof per node
            // for the purpose of all cases, nNodeDof[j] = 1, however left in for generality
        }

        // Local stiffness matrix (Ke) found by Gauss integration

        for (int f = 0; f < gaussorder; f++) {
            for (int g = 0; g < gaussorder; g++) {

                eta = GP[f];
                xi = GP[g];

                // below assignment specific to quadrilateral elements
                // manual assignment of N, although not required in calculations as GN is also manually assigned
//                N[0] = 0.25 * (1 - xi) * (1 - eta);
//                N[1] = 0.25 * (1 + xi) * (1 - eta);
//                N[2] = 0.25 * (1 + xi) * (1 + eta);
//                N[3] = 0.25 * (1 - xi) * (1 + eta);

                // manual assignment of GN
                GN[0] = -(1 - eta) * 0.25;
                GN[1] = -(1 - xi) * 0.25;
                GN[2] = (1 - eta) * 0.25;
                GN[3] = -(1 + xi) * 0.25;
                GN[4] = (1 + eta) * 0.25;
                GN[5] = (1 + xi) * 0.25;
                GN[6] = -(1 + eta) * 0.25;
                GN[7] = (1 - xi) * 0.25;

                // CBLAS matrix multiplication to obtain J
                // J = GN . eCoord
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2, 2, nnode_elem, 1.0, GN, 2, eCoord, 4, 0.0, J,
                            2);

                // manual calculation of determinate
                // May be replaced with lapack routine - however manual calculation is simple for quadrilateral elements
                // Lapack routine dgetrf produces LU decomposition of J: product of diagonals of U = det(A). Shown below
                // F77NAME(dgetrf)(2, 2, J, 2, ipiv, info);
                detJ = J[0] * J[3] - J[1] * J[2];

                // manual calculation of inverse
                // may also replace with lapack routine (getri), however inverse is trivial to calculate for 2x2
                invJ[0] = 1 / detJ * J[3];
                invJ[1] = -1 / detJ * J[2];
                invJ[2] = -1 / detJ * J[1];
                invJ[3] = 1 / detJ * J[0];

                // invJ . GN = B
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2, 4, 2, 1.0, invJ, 2, GN, 2, 0.0, B, 2);

                //-------- elemental stiffness matrix: Ke = Ke + B'*D*B*th*DetJ*Wi*Wj ------------
                // Ke = Ke + np.dot(np.dot(B.T, D), B) * th * DetJ * W[i] * W[j]

                // compute coefficient in front matrix multiplication: Ke = Ke + B'*D*B*k
                alpha = th * detJ * W[f] * W[g];

                // product inside brackets : B'. D = BtD
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 4, 2, 2, 1.0, B, 2, D, 2, 0.0, BtD, 4);

                // product outside brackets : k (BtD . B) + 1 * Ke = Ke
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 4, 2, alpha, BtD, 4, B, 2, 1.0, Ke, 4);


            }
        }

        // "inserting" the global element stiffness matrix into the global system stiffness matrix
        // loop though node DoF and insert into relevant row (f) and column (g)
        for (int f = 0; f < nnode_elem; f++) {
            for (int g = 0; g < nnode_elem; g++) {
                K[gDof[f] * nDof + gDof[g]] += Ke[f * nnode_elem + g];
            }
        }

    }


    // -----------------------------------------
    // **     NODAL FLUX & BOUNDARY CONDS    **
    // -----------------------------------------

    int fluxnodes[nFluxNodes];
    int tempnodes[nTempNodes];
    int nbe = nFluxNodes - 1;                   // number of flux elements
    double T[nDof];                             // nodal temperature vector
    double f[nDof];                             // nodal flux vector
    zeros_double(nDof, &(f[0]));                // zero overall flux and temperature vectors
    zeros_double(nDof, &(T[0]));
    double fq[2];                               // flux vector for nodes 1 & 2
    int node1;
    int node2;
    double n_bce[2];      // nodes 1 & 2 and flux value at edge
    float x1;
    float y1;
    float x2;
    float y2;     // nodal coordinates
    double leng;                                // element edge length
    double flux[1];                             // elemental flux as an array length 1 (for CBLAS)
    double N_1d[2];                             // shape function matrix
    float BC[nTempNodes * 2];                   // boundary condition nodes and corresponding temperatures
    int n_bc[4 * nbe];                          // elemental boundary condition matrix and node numbers
    int rDof = nDof - nTempNodes;

    // nodal locations of flux and constant-temperature vectors are case-specific
    switch (case_no) {

        // case 1
        case 1: {

            for (int i = 0; i < nFluxNodes; i++) {
                fluxnodes[i] = NodeTopo[i + (nnode_x - 1) * nnode_y];   // assign fluxnodes as nodes on right hand edge
            }

            for (int i = 0; i < nTempNodes; i++) {
                BC[i] = NodeTopo[i];                              // boundary condition nodes as left hand edge
                BC[i + nTempNodes] = T0;                          // boundary condition matrix for further manipulation
                tempnodes[i] = BC[i];
                T[NodeTopo[i]] = BC[i + nTempNodes];              // insert boundary temperature into global T vector
            }
            break;
        }

            // case 2: T_BOTTOM and Q_TOP
        case 2: {

            for (int i = 0; i < nFluxNodes; i++) {
                fluxnodes[i] = NodeTopo[i * nnode_y];             // assign fluxnodes as nodes on top edge
            }

            for (int i = 0; i < nTempNodes; i++) {
                BC[i] = NodeTopo[nnode_y-1 + i*nnode_y];            // boundary condition nodes as bottom edge
                BC[i + nTempNodes] = T0;
                tempnodes[i] = BC[i];
                T[NodeTopo[nnode_y-1 + i*nnode_y]] = BC[i + nTempNodes];   // insert boundary temperature into global T vector
            }
            break;
        }

            // case 3: T_LEFT and Q_BOTTOM
        case 3: {


            for (int i = 0; i < nFluxNodes; i++) {
                fluxnodes[i] = NodeTopo[(i + 1) * nnode_y - 1];       // assign fluxnodes as nodes on bottom edge
            }

            for (int i = 0; i < nTempNodes; i++) {
                BC[i] = NodeTopo[i];                               // boundary condition nodes as left edge
                BC[i + nTempNodes] = T0;
                tempnodes[i] = BC[i];
                T[NodeTopo[i]] = BC[i + nTempNodes];              // insert boundary temperature into global T vector
            }
            break;
        }
    }


    // flux for each element edge at the flux edge
    for (int i = 0; i < nbe; i++) {
        n_bc[i * 4] = fluxnodes[i];           // first row: Node 1
        n_bc[i * 4 + 1] = fluxnodes[i + 1];   // second row: Node 2
        n_bc[i * 4 + 2] = q;                  // third and fourth rows: flux in
        n_bc[i * 4 + 3] = q;
    }


    for (int i = 0; i<nbe; i++){

        // node number from boundary condition matrix, for each element edge
        node1 = n_bc[i*4];
        node2 = n_bc[i*4+1];
        n_bce[0] = n_bc[i*4+2];
        n_bce[1] = n_bc[i*4+3];

        // split nodal coordinates of element edge into separate variables
        x1 = Coord[node1];
        y1 = Coord[node1 + nnode];
        x2 = Coord[node2];
        y2 = Coord[node2 + nnode];

        leng = sqrt(    pow((x1 - x2), 2) + pow((y1 - y2), 2)   );      // element edge length

        detJ = leng/2;              // 1D Jacobian

        zeros_double(2, &(fq[0]));  // reset fq[] to zero for use in nodewise calculation

        // integrate in xi direction (1D)
        for(int j = 0; j < gaussorder; j++){

            xi = GP[j];     // 1D shape function point

            N_1d[0] = 0.5 * (1 - xi);   // populate shape function matrix
            N_1d[1] = 0.5 * (1 + xi);


            // flux = N_1d . n_bce
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 1, 1, 2, 1.0, N_1d, 1, n_bce, 2, 0.0, flux, 1);

            // nodal flux
            fq[0] = (fq[0] + W[j] * N_1d[0] * flux[0] * detJ * th);
            fq[1] = (fq[1] + W[j] * N_1d[1] * flux[0] * detJ * th);



        }

        f[node1] -= fq[0];      // subtract to make integral negative and insert into global flux vector
        f[node2] -= fq[1];

    }



        // ---------------------------
    // **  BOUNDARY CONDITIONS  **
    // ---------------------------

    // ----- Assembling global "Force" vector ------------
    // original Dof and temperature vector
//    int OrgDof[nDof];
//    zeros_int(nDof, &(OrgDof[0]));
//    int rDof = nDof;    // reduced number of Dof
//
//    for (int i = 0; i < nTempNodes; i++){
//
//        OrgDof[tempnodes[i]] = - 1;
//        T[tempnodes[i]] = BC[i + nTempNodes];
//
//    }
//    // update reduced degrees of freedom
//    rDof -= nTempNodes;
//
//    int RedDof[rDof];
//    int counter1 = 0;
//    zeros_int(rDof, &(RedDof[0]));
//
//    for (int i = 0; i < nDof; i++){
//
//        if (OrgDof[i] == 0){
//            OrgDof[i] = counter1;
//            RedDof[counter1] = i;
//            counter1++;
//        }
//    }




    // ---------------------------
    // **  PARTITION MATRICES  **
    // ---------------------------

    const int E_n = nTempNodes;            // size of mask_E returning true values, used later for size of matrix partitions
    const int Ef_n = nDof - nTempNodes;    // size of mask_E returning false values
    int mask_E[nnode];                 //  mask used for partitioning
    int mask_Ef[nnode];              //  elementwise inverse is also required
    double T_E[E_n];                 // temperature vector of boundary coundition nodes
    double f_F[Ef_n];                // flux vector of flux input nodes
    double T_F[Ef_n];
    int T_Ecount = 0;
    int f_Fcount = 0;
    double K_EE[E_n * E_n];          // partition of K corresponding to fixed temperature nodes
    double K_FF[Ef_n * Ef_n];        // partition of K corresponding to flux nodes
    double K_EF[E_n * Ef_n];
    double rhs[Ef_n];                // rhs of system to be solved for
    int count_E = 0;
    int count_Ef = 0;
    zeros_int(nnode, &(mask_E[0]));

    for(int i = 0; i < nnode; i++){
        mask_Ef[i] = 1;
    }
    //  create boolean mask  -
//  ---------> NOTE: this mask is correct for all cases, however the result after is incorrect
//    for(int i = 0; i < nTempNodes; i++){
//    zeros_int(nnode, &(mask_E[0]));
//
//    for(int i = 0; i < nnode; i++){
//        mask_Ef[i] = 1;
//    }
//
//        mask_E[tempnodes[i]] = 1;
//        mask_Ef[tempnodes[i]] = 0;
//    }
//    int countKEE = 0;
//    int countKEF = 0;
//    int countKFF = 0;
//    for (int i = 0; i < nnode; i++) {
//        for (int j = 0; j < nnode; j++) {
//
//            if (mask_E[i] == 1 && mask_E[j] == 1){
//                K_EE[countKEE] = K[nnode * i + j];
//                countKEE++;
//            }
//            else if (mask_E[i] == 1 && mask_E[j] != 1){
//                K_EF[countKEF] = K[nnode * i + j];
//                countKEF++;
//            }
//            else if (mask_E[i] != 1 && mask_E[j] != 1){
//                K_FF[countKFF] = K[nnode * i + j];
//                countKFF++;
//            }
//        }
//    }



    // this mask works for case 1 only
    for(int i = 0; i < nDof; i++) {
        mask_E[i] = false;                       // initialise each mask entry as false
        mask_Ef[i] = true;
        for (int j = 0; j < nTempNodes; j++) {
            if (i == tempnodes[j]) {
                mask_E[i] = true;                // loop through TempNodes and return true if any value 1:nDof matches
                mask_Ef[i] = false;              // invert for mask_Ef
            }
        }
    }

    // ----- populating partitioned matrices ------------------
    //set size based on mask_E. Loop though partition and insert into relevant row (j) and column (i) from K
    //since global stiffness may be split into consecutive 'chunks'

    // K_EE -- python: K_EE = K[np.ix_(mask_E, mask_E)]
    for (int i = 0; i < E_n; i++) {
        for (int j = 0; j < E_n; j++) {
            K_EE[i * E_n + j] =  K[i * nDof + j];
        }
    }

    // K_FF = K[np.ix_(~mask_E, ~mask_E)]
    // K_EF = K[np.ix_(mask_E, ~mask_E)]
    for (int i = 0; i < Ef_n; i++) {

        for (int j = 0; j < Ef_n; j++) {
            K_FF[i * Ef_n + j] =  K[(i+E_n) * nDof + (j+E_n)];
        }

        for (int j = 0; j < E_n; j++) {
            K_EF[i * E_n + j] =  K[(i+E_n) * nDof + (j)];
        }
    }


    for(int i = 0; i < nDof; i++) {

        if (mask_E[i] == 1) {                         // for true values of mask (where T[i] == T0, populate T_E
            T_E[T_Ecount] = T[i];
            T_Ecount++;
        } else {                                  // for false values, populate f_F
            f_F[f_Fcount] = f[i];
            f_Fcount++;
        }
    }


    // ---------------------------
    // **  CONSTRUCT SOLUTION  **
    // ---------------------------

    // rhs = f_F - K_EF' . T_E
    cblas_dgemv(CblasColMajor, CblasTrans, E_n, Ef_n, 1.0, K_EF, E_n, T_E, 1, 0.0, rhs, 1);
    for (int i = 0; i < Ef_n; i++) {
       rhs[i] =  f_F[i] - rhs[i];
    }


    int ipiv[Ef_n];     // pivot points of solution
    int info=0;                     // lapack routine - provides basic information whether routine is correctly executed.


    F77NAME(dgesv)(Ef_n, 1, K_FF, Ef_n, ipiv, rhs, Ef_n, info);        // solve system of linear equations
    cblas_dcopy(Ef_n, rhs, 1, T_F, 1);                                     // copy rhs to T_F


    for (int i = 0 ; i < Ef_n ; i++) {
        T[i + E_n] = T_F[i];
    }


    double f_E[nTempNodes];
    // calculate dot products np.dot(K_EE, T_E) + np.dot(K_EF, T_F). Sum result of first one with second using beta = 1
    cblas_dgemv(CblasColMajor, CblasNoTrans, E_n, E_n, 1.0, K_EE, E_n, T_E, 1, 0.0, f_E, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, E_n, Ef_n, 1.0, K_EF, E_n, T_F, 1, 1.0, f_E, 1);


    for (int i = 0; i < E_n; i++){
        f[i] = f_E[i];
    }

    for (int i = 0; i < Ef_n; i++){
        f[i + E_n] = f_F[i];
    }


    MPI_Finalize();

    // call vtkwrite function to write to output file
    vtkwrite(nnode, &(Coord[0]), nelem, nnode_elem, &(ElemNode[0]), &(T[0]));


    return 0;


}



