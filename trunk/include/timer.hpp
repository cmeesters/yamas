/***************************************************
 *                                                 *
 * YAMAS - Yet Another Meta-Analysis               *
 *                                                 *
 * This file contains the YAMAS timing facilities. *
 * To be used for internal development purposes,   *
 * only.                                           *
 * This file is not part of the yamas namespace    *
 * and intended for general purpose profiling.     *
 *                                                 *
 **************************************************/

// $Rev: 922 $
// $LastChangedDate: 2012-08-01 12:52:00 +0200 (Mi, 01 Aug 2012) $
// $Author: meesters $

#include <ctime>
#include <iostream>
using namespace std;

#ifndef TIMERHPP
#define TIMERHPP

class Timer {
    public:
        Timer(int number = 0) : n(number), 
                                start(clock()) {}
        ~Timer();
    private:
        int n;
        clock_t start;
};
 
inline Timer::~Timer() {
    cout << "timer " << n << " took "
         << double(clock() - start) / CLOCKS_PER_SEC
         << " seconds" << endl;
}
#endif //TIMER

/* The following part shows a brief, but complete usage
   example:

int main() {
   vector<int> vi;
   {
       Timer t0(0);
       for (int i = 0; i < 50000; ++i) {
           vi.insert(vi.begin(), rand());
       }
   } // implicit call of ~Timer
   {
       Timer t1(1);
       vector<int>::const_iterator imax = max_element(vi.begin(), vi.end());
       cout << "The largest element is " << *imax << endl;
       vector<int>::const_iterator imin = min_element(vi.begin(), vi.end());
       cout << "The smallest element is " << *imin << endl;
   } // implicit call of ~Timer
}
*/
