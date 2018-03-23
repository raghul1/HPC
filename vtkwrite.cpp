#include <fstream>
#include <iomanip>
using namespace std;

void vtkwrite(int nnode, double* Coord, int nelem, int nnode_elem, int* ElemNode, double* T){

    ///@brief write in ASCII to a  .vtk file, which may be used in paraview
    ///@brief syntax followed using python ASCII output

    ///@param nnode         total number of nodes
    ///@param Coord         coordinates (x,y) of each node
    ///@param nelem         total number of elements
    ///@param nnode_elem    number of elements per node- for quadrilateral elements this is restricted to 4
    ///@param ElemNode      5 columns: 1x element number, 4x corners for a quadrilateral element
    ///@param T             array of temperature values corresponding to nodal temperatures

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
