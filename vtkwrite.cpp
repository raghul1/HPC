#include <fstream>
#include <iomanip>
using namespace std;

void vtkwrite(int nnode, double* Coord, int nelem, int nnode_elem, int* ElemNode, double* T){

    // ---------------------------
    // **    VTK FILE WRITE    **
    // ---------------------------

    // write to ASCII .vtk file, which may be used in paraview
    // syntax followed using python ASCII output

    ofstream vOut("output.vtk", ios::out | ios::trunc);      // output to file, overwrite if already present

    // 'preamble'
    vOut << "# vtk DataFile Version 4.0" << endl;
    vOut << "vtk output" << endl;
    vOut << "ASCII" << endl;
    vOut << "DATASET UNSTRUCTURED_GRID" << endl;

    // write nodal coordinates as one line. Coordinates separated by single spaces and each node is separated by 0.0
    vOut << "POINTS " << nnode << " " << "double" << endl;
    for(int i = 0; i < nnode; i++){
        vOut << Coord[i] << " " << Coord[i + nnode] << " " << 0.0 << " ";
    }
    vOut << endl;

    // write nodal numbers for each element in form (total nodes per element, upper left, upper right, lower right, lower left)
    vOut << "CELLS " << nelem << " " << nelem * 5 << endl;
    for(int i = 0; i < nelem; i++) {
        vOut  << nnode_elem << " " << ElemNode[i + nelem] << " " << ElemNode[i + nelem * 2] << " " << ElemNode[i + nelem * 3] << " "
              << ElemNode[i + nelem * 4] << endl;
    }

    // cell types
    vOut << "CELL_TYPES " << nelem << endl;
    for(int i = 0; i < nelem; i++) {
        vOut << 9 << endl;
    }

    // output formatting - temperature
    vOut << "POINT_DATA " << nnode << endl;
    vOut << "FIELD FieldData " << 1 << endl;
    vOut << "disp " << 1 << " " << nnode << " " << "double" << endl;
    for(int i = 0; i < nnode; i++) {
        vOut << T[i] << " ";
    }
    vOut << endl;
}
