#include <math.h>
#include <array>
#include <deque>
#include <iostream>
#include <limits.h>
#include "Neighbor.h"

using namespace std;

// Constructors
Neighbor::Neighbor() : GrowArray() {};

Neighbor::Neighbor(int i) : GrowArray(i) {};

Neighbor::Neighbor(int i, int j) : GrowArray(i, j) {};

// Destructor
Neighbor::~Neighbor() {};

// Generator
void Neighbor::GeneratePeriodic() {
    // If it's a left or right shift
    if ((lrshift != 0) && (udshift == 0)) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < array[i].size(); ++j) {
                if (j + lrshift < 0) {
                    array[i][j] = array[i].size() + j + lrshift;
                } else {
                    array[i][j] = (j + lrshift) % array[i].size();
                }
            }
        }
    }
    else if ((lrshift == 0) && (udshift != 0)) {
        if (array.size() != 1){
            for (int i = 0; i < rows; ++ i) {
                for (int j = 0; j < array[i].size(); ++j) {
                    int idx2test = i + udshift;
                    if (checkExist(idx2test, j) == true) { //if it's a valid array point
                        array[i][j] = idx2test;
                    } else {
                        while (checkExist(idx2test, j) == false) {
                            if (idx2test < 0) {  // circles back to bottom row of the array
                                idx2test = rows - 1;
                            }
                            else if (idx2test > rows - 1) {  // circles back to 0th row of array
                                idx2test = idx2test % rows;
                            } else {
                                idx2test = idx2test + udshift;
                            }
                        }
                        array[i][j] = idx2test;
                    }
                }
            }
        }
        else {
             for (int i = 0; i < rows; ++ i) {
                for (int j = 0; j < cols; ++j) {
                    array[i][j] = 0;
                }
            }
        }
    }
    else if ((lrshift == 0) && (udshift == 0)) {
        cout << "Need to define a shift for Neighbors" << endl;
    } else {
        cout << "Error: Only one directional shift supported at this time" << endl;
    }
}

//Growth functions  //Note: extend shouldn't do anything, just there as placeholder to match form of parent function (necessary to override)
void Neighbor::grow1D(bool extend) { 

    //Add one space to each row
    int numRows = array.size();
    for (int i =0; i < numRows; ++i){
        array[i].push_back(0);
    }
    cols = cols+1;
    
    //Re-generate new values
    GeneratePeriodic();
}

void Neighbor::growTrap(bool vertextend, bool horizextend) {
    //First, add a column to the left
    for (unsigned short i = 0; i < rows; ++i) {
        array[i].push_front(0);
    }

    //Then, add new rows
    deque<unsigned short> newDeque(1, 0);
    array.push_front(newDeque);
    array.push_back(newDeque);

    //Update indices
    rows = rows + 2;
    cols = cols + 1;

    GeneratePeriodic();

    if (rows != array.size()) {
        cout << "Error in row indexing" << endl;
    }
}

