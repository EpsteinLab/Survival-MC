#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <array>
#include <ctime>
#include <time.h>
#include <random>
#include <string>
#include <fstream>
#include <limits.h>
#include <algorithm>

#include "include/GrowArray.h"
#include "include/ZArray.h"
#include "include/Irid.h"
#include "include/Neighbor.h"
#include "include/HitTracker.h"

using namespace std;

int main()
{
///////////////////////////////////////////////////////////////////////////////////////////////
    // Please change the path to where you want to store simulations
    // This MUST be done to get output .csv files
    std::string saveSim = "file/path/to/Simulation/Folder;

////////////////////////////////////////////////////////////////////////////////////////////////

    // CHANGE SIMULATION AND "REACTION" PARAMETERS HERE:

    // Initial lattice dimensions
    short unsigned int r = 400; //rows
    short unsigned int c = 1; //columns
    

    // Physical and "Reaction" Parameters
    short unsigned int h = 30; //Long-range interaction distance
    float bx = 1; //natural birth rate of xanthophores
    float bm = 7; //natural birth rate of melanophores
    float dx = 0; //natural death rate of xanthophores
    float sm = 1; //killing of xanthophores by melanophores
    float sx = 1; //killing of melanophores by xanthophores
    float dm = 9; //death of melanophores with non-xanthophores h distance away
    float dmx = 0; //death rate of melanophores with xanthophores h distance away

    // Simulation Parameters
    unsigned long long int stepsPerGrowth = 1e7; //Number of MC iterations before growth. We recommend 10^7 as a standard.
                // There is no noticeable change if the number is increased above 10^7, though the simulations take much longer.
                // Below 10^7 the simulations begin to get messy/don't settle into stable patterns
    unsigned long long int growSteps = stepsPerGrowth * r; //Number of MC steps - these settings grow the lattice to a square
    unsigned long long int imgPerSim = growSteps / 1000; //Number of .csv simulation files exported
    unsigned long long int endSteps = 0 * stepsPerGrowth; //Allows for simulation to carry on past the end of the growth
    unsigned long long int beginSteps = 0 * stepsPerGrowth; //Allows for simulation to iterate before growth begins


    //Ablation Parameters  
    int ablateWhen1000 = 0; //Change this value for when you want the simulation to ablate
                // If you do not want ablation to occur, set this to zero
                // If you do want an ablation, change this value to 1000*the fraction of growth you want to ablate at
                // For example, if you want the system to ablate halfway through growth, set ablateWhen1000 = 500;
    float widthAblate = 0.75; //These tell us to ablate this proportion of the rows widths
    float heightAblate = 0.75; //This tells us wot ablate this proportion of the columns height


    // Create Initial Arrays
    ZArray zebra(r, c); // The main array where it stores the values of the chromatophores
    zebra.BlankArray(); // Can change this function to a different function for initialization, though it has no noticeable long-term effect

    Irid ir(r, c); // Iridophore array, for guiding the patterns
    ir.Blank(); // Can change this function to add iridophores into the simulation
                //Iridophores make it so melanophores cannot form at that lattice node



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This section is necessary to run the simulations, but should not be changed
    // Instead, use the parameters above to change what is happening here

    short unsigned int r0 = r;  //Initial Condtions (to export - not used in simulations)
    short unsigned int c0 = c;
    
    // Calculate probabilities
    float rate_sum = bx + bm + dx + dmx + dm + sm + sx;
    float bx_p = bx / rate_sum;
    float bm_p = bm / rate_sum;
    float dx_p = dx / rate_sum;
    float dmx_p = dmx / rate_sum;
    float sm_p = sm / rate_sum;
    float sx_p = sx / rate_sum;
    float dm_p = dm / rate_sum;

    // Calculates total number of steps
    unsigned long long int fullSteps = endSteps + beginSteps + growSteps;

    // Tells system when to ablate
    list<unsigned long long int> ablateSpots;
    if (ablateWhen1000 != 0) {
        ablateSpots.push_back(((growSteps * ablateWhen1000) / 1000) + beginSteps);
    } else {
        ablateSpots.push_back(0);
    }

    // Create Neighbor Arrays
    Neighbor up(r, c); // A neighbor matrix to pull array values from
    up.udshift = -1;
    up.GeneratePeriodic();

    Neighbor down(r, c);
    down.udshift = 1;
    down.GeneratePeriodic();

    Neighbor left(r, c);
    left.lrshift = -1;
    left.GeneratePeriodic();

    Neighbor right(r, c);
    right.lrshift = 1;
    right.GeneratePeriodic();

    HitTracker hitCnt;

//////////////////////////////////////////////////////////////////////////////////////////////////////

    // MONTE CARLO SIMULATION ALGORITHM

    //Run Simulation by looping through each event. 
    mt19937 generator ((int) time(0));
    uniform_real_distribution<float> dis(0.0, 1.0);
    srand(time(0));

    unsigned long long hitCounter = 0;
    for (unsigned long long int iter = 0; iter < fullSteps; ++iter) {
        c = zebra.getCols();
        r = zebra.getRows();
        float proc = dis(generator); // Chooses random process
        short int i = rand() % r; // Chooses random lattice point
        short int j = rand() % c;
        
        if (zebra.checkExist(i, j) == true) {
            if (proc < dm_p) { //death of melanophores with non xanthophores h distance away
                if (zebra.array[i][j] == 2) {
                    float angle = dis(generator) * 2 * M_PI;
                    float cosangle = cos(angle);
                    float sinangle = sin(angle);
                    short int inew = rint(sinangle * h + i);
                    short int jnew = rint(cosangle * h + j);
                    if (zebra.checkExist(inew, jnew) == false) {
                        // First, check and make sure indicies are in area.
                        while (inew < 0) {
                            inew = zebra.getRows() + inew;
                        } 
                        if (inew >= zebra.getRows()) {
                            inew = inew % zebra.getRows();
                        }
                        unsigned numCol = zebra.array[inew].size();
                        if (jnew < 0) {
                            jnew = (numCol - 1 + jnew) % numCol;
                        }
                        else if (jnew >= numCol) {
                            jnew = jnew - numCol;
                        }
                    }
                    if ((zebra.array[inew][jnew] == 0) || (zebra.array[inew][jnew] == 2)) {
                        zebra.array[i][j] = 0;  // Die out
                        hitCounter = hitCounter + 1;
                    }
                    
                }
            }

            else if ((proc < dm_p + dmx_p) && (proc >= dm_p)) { //death of melanophore with xanthophore h away
                if (zebra.array[i][j] == 2) {
                    float angle = dis(generator) * 2 * M_PI;
                    float cosangle = cos(angle);
                    float sinangle = sin(angle);
                    short int inew = rint(sinangle * h + i);
                    short int jnew = rint(cosangle * h + j);
                    if (zebra.checkExist(inew, jnew) == false) {
                        // First, check and make sure indicies are in area.
                        while (inew < 0) {
                            inew = zebra.getRows() + inew;
                        } 
                        if (inew >= zebra.getRows()) {
                            inew = inew % zebra.getRows();
                        }
                        unsigned numCol = zebra.array[inew].size();
                        if (jnew < 0) {
                            jnew = (numCol - 1 + jnew) % numCol;
                        }
                        else if (jnew >= numCol) {
                            jnew = jnew - numCol;
                        }
                    }
                    if (zebra.array[inew][jnew] ==1)  {
                        zebra.array[i][j] = 0;  // Die out
                        hitCounter = hitCounter + 1;
                    }
                    
                }
            }

            else if ((proc < dm_p + dmx_p + sm_p) && (proc >= dm_p + dmx_p)) { //Tests to see if local melanophore will kill xanthophore
                if (zebra.array[i][j] == 1) {
                    switch(rand() % 4) {
                        case 0: { //Above Neighbor
                            if ((up.array[i][j] != SHRT_MAX) && (up.checkExist(i, j) == true)) {
                                if (zebra.array[up.array[i][j]][j] == 2) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        case 1: { //Below Neighbor
                            if ((down.array[i][j] != SHRT_MAX) && (down.checkExist(i, j) == true)) {
                                if (zebra.array[down.array[i][j]][j] == 2) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;    
                            }
                        }
                        case 2: { //Left Neighbor
                            if ((left.array[i][j] != SHRT_MAX) && (left.checkExist(i, j) == true)) {
                                if (zebra.array[i][left.array[i][j]] == 2) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        case 3: { //Right Neighbor
                            if ((right.array[i][j] != SHRT_MAX) && (right.checkExist(i, j) == true)) {
                                if (zebra.array[i][right.array[i][j]] == 2) {
                                    zebra.array[i][j] = 0;  
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        default: {
                            break;
                        }
                    }
                }
            } 

            else if ((proc < dm_p + dmx_p + sm_p + sx_p) && (proc >= dm_p + dmx_p + sm_p)) { //Tests to see if local xanthophore will kill melanophore
                if (zebra.array[i][j] == 2) {
                    switch(rand() % 4) {
                        case 0: { //Above Neighbor
                            if ((up.array[i][j] != SHRT_MAX) && (up.checkExist(i, j) == true)) {
                                if (zebra.array[up.array[i][j]][j] == 1) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        case 1: { //Below Neighbor
                            if ((down.array[i][j] != SHRT_MAX) && (down.checkExist(i, j) == true)) {
                                if (zebra.array[down.array[i][j]][j] == 1) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;    
                            }
                        }
                        case 2: { //Left Neighbor
                            if ((left.array[i][j] != SHRT_MAX) && (left.checkExist(i, j) == true)) {
                                if (zebra.array[i][left.array[i][j]] == 1) {
                                    zebra.array[i][j] = 0;
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        case 3: { //Right Neighbor
                            if ((right.array[i][j] != SHRT_MAX) && (right.checkExist(i, j) == true)) {
                                if (zebra.array[i][right.array[i][j]] == 1) {
                                    zebra.array[i][j] = 0;  
                                    hitCounter = hitCounter + 1;
                                }
                                break;
                            }
                        }
                        default: {
                            break;
                        }
                    }
                }
            }

            else if ((proc < dm_p + dmx_p + sm_p + sx_p + bx_p) && (proc >= dm_p + dmx_p + sm_p + sx_p)) { //Tests for birth of xanthophore
                if (zebra.array[i][j] == 0) {
                    zebra.array[i][j] = 1;
                    hitCounter = hitCounter + 1;
                }
            }

            else if ((proc < dm_p + dmx_p + sm_p + sx_p + bx_p + bm_p) && (proc >= dm_p + dmx_p + sm_p + sx_p + bx_p)) { //Tests for birth of melanophore
                if ((zebra.array[i][j] == 0) && (ir.array[i][j] == 0)) {
                    zebra.array[i][j] = 2;
                    hitCounter = hitCounter + 1;
                }
            }

            else if ((proc < dm_p + dmx_p + sm_p + sx_p + bx_p + bm_p + dx_p) && (proc >= dm_p + dmx_p + sm_p + sx_p + bx_p + bm_p)) { //Tests for death of xanthophore
                if (zebra.array[i][j] == 1) {
                    zebra.array[i][j] = 0;
                    hitCounter = hitCounter + 1;
                }
            }

            else {
                std::cout << "Error: Incorrect number generated" << endl;
            }
        }

        
        //Check to see if the system needs to be ablated
        if (find(ablateSpots.begin(), ablateSpots.end(), iter) != ablateSpots.end()) {
            zebra.Ablate(widthAblate, heightAblate);
            unsigned long itnum = floor(iter / imgPerSim);
            string iter_num = to_string(itnum);
            if (iter_num.length() == 1){
                iter_num = "0000" + iter_num;
            } else if (iter_num.length() == 2){
                iter_num = "000" + iter_num;
            } else if (iter_num.length() == 3){
                iter_num = "00" + iter_num;
            }else if (iter_num.length() == 4){
                iter_num = "0" + iter_num;
            }

            string outname = saveSim + "/img_" + iter_num + "_ablation.csv";
            zebra.Rectangular_Export(outname);

            hitCnt.hitCounts.push_back(hitCounter);
            hitCnt.iterValues.push_back(iter);
        }
        

        // Export image regularly
        if (iter % imgPerSim == 0) {
            string iter_num = to_string(iter / imgPerSim);
            if (iter_num.length() == 1){
                iter_num = "0000" + iter_num;
            } else if (iter_num.length() == 2){
                iter_num = "000" + iter_num;
            } else if (iter_num.length() == 3){
                iter_num = "00" + iter_num;
            }else if (iter_num.length() == 4){
                iter_num = "0" + iter_num;
            }

            string outname = saveSim + "/img_" + iter_num + ".csv";
            zebra.Rectangular_Export(outname);

            hitCnt.hitCounts.push_back(hitCounter);
            hitCnt.iterValues.push_back(iter);
        }

        // Perform Growth
        if ((iter % stepsPerGrowth == 0) && (iter <= growSteps + beginSteps) && (iter >= beginSteps)) {
            zebra.grow1D(false);
            ir.grow1D(true);  //This should ALWAYS be true: the iridophores guide the pattern, so need to extend
            up.grow1D(false);
            down.grow1D(false);
            left.grow1D(false);
            right.grow1D(false);
        }
        
    }   
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    // EXPORT A CSV OF THE CONDITIONS
    
    // Check if there are iridophores
    string iridAns;
    if (ir.checkExist()) {
        iridAns = "True";
    } else {
        iridAns = "False";
    }

    // Export conditions to csv
    string csvCondTitle = saveSim + "/conditions.csv";
    ofstream csvfile;
    csvfile.open(csvCondTitle);

    csvfile << "Initial_Rows" << "," << "Initial_Columns" << "," << "h" << ","; 
    csvfile << "Growing_Steps" << "," << "Steps_per_Growth" << ",";
    csvfile << "Images" << "," << "Num_Hits" << "," ;
    csvfile << "Final_Rows" << "," << "Final_Columns" << "," << "Irid_Exist" << ",";
    csvfile << "bx" << "," << "bm" << "," << "dx" << "," << "dm" << "," << "dmx" << "," ;
    csvfile << "sm" << "," << "sx" << "," << "Boundary Conditions" << ","<< "Ablated"<< endl;

    csvfile << to_string(r0) << "," << to_string(c0) << "," << to_string(h) << ",";
    csvfile << to_string(growSteps) << "," << to_string(stepsPerGrowth) << ",";
    csvfile << to_string(imgPerSim)<< "," << to_string(hitCounter) << ",";
    csvfile << to_string(r) << "," << to_string(c) << "," << iridAns << "," ;
    csvfile << to_string(bx) << "," << to_string(bm) << "," << to_string(dx) << "," << to_string(dm) << "," << to_string(dmx) << "," ;
    csvfile << to_string(sm) << "," << to_string(sx) << "," << "Periodic" <<",";
    
    for (auto const &v : ablateSpots) {
        csvfile << to_string(v)<< ";";
    }


    csvfile << endl;
    

    csvfile.close();

    // Export hit table to csv
    std::string hitName = saveSim + "/hitNums.csv";
    hitCnt.export2csv(hitName);


   std::cout << "Completed Simulation" << endl;
}
