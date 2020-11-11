#include <math.h>
#include <fstream>
#include <cstring>
#include <deque>

#include "GrowArray.h"

using namespace std;

// Constructors
GrowArray::GrowArray() : rows(0), cols(0) {};

GrowArray::GrowArray(int i) : rows(i), cols(i) {  //with one argument, makes square
    deque<unsigned short> blankRow;
    for (unsigned short k = 0; k < cols; ++k){  //Constructs blank row deque
        blankRow.push_back(0);
    }
    for (unsigned short ll = 0; ll < rows; ++ll){  //fills array with blank deque
        array.push_back(blankRow);
    }
}

GrowArray::GrowArray(int i, int j) : rows(i), cols(j) {
    deque<unsigned short> blankRow;
    for (unsigned short k = 0; k < cols; ++k){  //Constructs blank row deque
        blankRow.push_back(0);
    }
    for (unsigned short ll = 0; ll < rows; ++ll){  //fills array with blank deque
        array.push_back(blankRow);
    }
}

// Destructor
GrowArray::~GrowArray() {}

// Growth Functions
void GrowArray::grow1D(bool extend) { // will grow cols by 1
    
    //Add one more to each column deque
    int numRows = array.size();
    for (int i =0; i < numRows; ++i){
        if (extend == false){
            array[i].push_back(0);
        } else {
            unsigned short oldVal = array[i].back();
            array[i].push_back(oldVal);
        }
    }
    //Increases column size
    cols = cols + 1;
}


void GrowArray::growTrap(bool vertextend, bool horizextend) {
    //Grows the system by a singular column with an extra row and column at the end, so it expands the size

    //First, add new column to the existing rows
    for (unsigned short i = 0; i < rows; ++i) {
        if (horizextend == false) {
            array[i].push_front(0);
        } else {
            array[i].push_front(array[i].front());
        }
    }

    //Then, add a new deque above and below
    if (vertextend == false) {
        deque<unsigned short> newDeque(1, 0);
        array.push_front(newDeque);
        array.push_back(newDeque);
    } else {
        deque<unsigned short> topDeque(1, array[0].front());
        array.push_front(topDeque);
        deque<unsigned short> bottomDeque(1, array[rows-1].front());
        array.push_back(bottomDeque);
    }

    //Update indices
    rows = rows + 2;
    cols = cols + 1;

    if (rows != array.size()) {
        throw invalid_argument("Error in row indexing");
    }
}


// Access Protected Dimensions
bool GrowArray::checkExist(int i, int j) {
    if ((i < rows) && (i >= 0)) {
        if ((j >= 0) && (j < array[i].size())) {
            return true;
        } else {
            return false;
        }
    } else {
    return false;
    }
}

unsigned short int GrowArray::getCols() {
    return cols;
}

unsigned short int GrowArray::getRows() {
    return rows;
}

// Cleanup and Export
void GrowArray::export2csv(string name) {
    ofstream tmpfilee;
    tmpfilee.open(name);
    for (int i=0; i<rows; ++i) {
        tmpfilee << to_string(array[i][0]);
        for (int j=1; j<cols; ++j) {
            tmpfilee << "," << to_string(array[i][j]);    
        }
        tmpfilee << endl;
    }
    tmpfilee.close();
}

void GrowArray::Rectangular_Export(string name) {
    //Rectangularize
    deque<deque<unsigned short>> output;
    for (unsigned short i = 0; i < rows; ++i){
        deque<unsigned short> thisRow = array[i];
        for (unsigned short j = thisRow.size() -1; j < cols; ++j) {
            thisRow.push_back(100);
        }
        output.push_back(thisRow);
    }

    //Export
    ofstream tmpfilee;
    tmpfilee.open(name);
    for (int i=0; i<rows; ++i) {
        tmpfilee << to_string(output[i][0]);
        for (int j=1; j<cols; ++j) {
            tmpfilee << "," << to_string(output[i][j]);    
        }
        tmpfilee << endl;
    }
    tmpfilee.close();
}



