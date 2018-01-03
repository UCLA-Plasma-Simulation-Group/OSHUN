/*!\brief  Functors for various time-integration methodds - Declarations
 * \author PICKSC
 * \file   clock.h
 *
 * Includes functors for fully explicit, implicit B, implicit E
 *
 */
#ifndef OSHUN_CLOCK_H
#define OSHUN_CLOCK_H

//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can
//  be incremented;

class Stepper {
public:
//      Constructor
    Stepper(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails);
    ~Stepper();

    double check_temperature(const State1D& Ystar, const State1D& Y);
    double check_flds(const State1D& Ystar, const State1D& Y, double& maxval);
    double check_flds(const State2D& Ystar, const State2D& Y, double& maxval);
    double check_js(const State1D& Ystar, const State1D& Y);

    bool update_dt(const State1D& Y_old, const State1D& Ystar, State1D& Y_new);
    bool update_dt(const State2D& Y_old, const State2D& Ystar, State2D& Y_new);

    Stepper& operator++();

    double dt() {return _dt;}
    double nextdt() {return dt_next;}
    double time() {return current_time;}

private:

    double current_time, dt_next, _dt;
    double atol, rtol, acceptability, err_val;
    
    size_t failed_steps, max_failures;

    int overall_check;

    size_t Nbc;
    int world_rank, world_size;

    double* acceptabilitylist;
};
//--------------------------------------------------------------

//
//**************************************************************
#endif //OSHUN_CLOCK_H

