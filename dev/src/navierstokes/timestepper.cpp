#include "field.h"
#include "timestepper.h"

TimeStepper::TimeStepper( double Dt, double T, double cfl )
{
    curtime = dt = 0; maxdt = Dt; tstop = T; nstep = 0; CFL = cfl;
    if (CFL>1 || CFL<=0) {
        throw FieldErr("Operator::ConvectionTimeStepper::ConvectionTimeStepper",
                       "CFL out of range");
    }
    history.push_back(0);
}

double TimeStepper::operator ++ (int)
{
    // save old time;
    double old_time = curtime;
    // update time, step counter and history
    curtime += dt;
    nstep ++;
    if (nstep < 0) {
        throw FieldErr("TimeStepper::operator++","nstep should not be negative");
    }
    // compute next timestep
    if ((unsigned int)nstep >= history.size()) {
        history.push_back(curtime);
        // compute next timestep
        dt = CFL * NextTimestep();
        // consider timestep limit and stop time
        if (dt > maxdt) dt = maxdt;
        if (curtime + dt >= tstop - 1E-10) {
            dt = tstop - curtime + 1E-10;
        }
    }
    else dt = history[nstep+1] - curtime;

    return old_time;
}

double TimeStepper::operator -- (int)
{
    if (nstep <= 0) {
        throw FieldErr("TimeStepper::operator --","Big Bang Reached");
    }
    // save old time
    double old_time = curtime;
    // update step count, time and timestep
    nstep --;
    curtime = history[nstep];
    dt   = history[nstep+1] - curtime;

    return old_time;
}

double TimeStepper::NextTimestep(void) {
    return maxdt / CFL;
}

double TimeStepper::SetStep(unsigned int n)
{
    if (n < 0 || n >= history.size()) {
        throw FieldErr("TimeStepper::SetStep","Out of time range");
    }
    nstep = n;
    curtime = history[nstep];
    if (n == history.size()-1) {
        // recompute next timestep
        dt = CFL * NextTimestep();
        // consider timestep limit and stop time
        if (dt > maxdt) dt = maxdt;
        if (curtime + dt >= tstop - 1E-10) {
            dt = tstop - curtime + 1E-10;
        }
    }
    else dt = history[nstep+1] - curtime;

    return curtime;
}

