#ifndef __NEIGHBOR_H_INCLUDED__
#define __NEIGHBOR_H_INCLUDED__

#include "GrowArray.h"

class Neighbor : public GrowArray {

public:
    // Constructors and Destructor
    Neighbor();
    Neighbor(int i);
    Neighbor(int i, int j);
    ~Neighbor();

    // Generators
    void GeneratePeriodic();

    // Growth - Need to override from parent class since we use these differently
    void grow1D(bool extend);  // Should be the same essentially as grow1ColBack
    void growTrap(bool vertextend, bool horizextend);


    // Shifts
    short lrshift {0};
    short udshift {0};

};


#endif