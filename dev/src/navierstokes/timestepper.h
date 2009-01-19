#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <vector>

using namespace std;

/* ------------------------------------------------------------------ */
/*                          Basic Time Stepper                        */
/* ------------------------------------------------------------------ */

class TimeStepper {
    private:
    double curtime, maxdt, tstop, dt, CFL;
    int nstep;
    vector<double> history;

    virtual double NextTimestep(void);

    public:
    TimeStepper( double dt, double T, double CFL );
    virtual ~TimeStepper() {}
    operator double& ( ) { return curtime;    }
    double TimeAt(int i) { return history[i]; }
    double Dt(void)      { return dt;         }
    int NStep(void)      { return nstep;      }
    void   SetStopTime(double T) { tstop = T; }
    double operator ++ (int);
    double operator -- (int);
    double SetStep(unsigned int n);
};

/* ------------------------------------------------------------------ */
/*                        Solution History class                      */
/* ------------------------------------------------------------------ */

template <class SOLN>
class SolnHistory {
    private:
    vector<SOLN> memory;
    vector<int>  tsteps;

    public:
};

#endif
