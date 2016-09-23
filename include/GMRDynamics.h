/*
 * GMRDynamics.h
 *
 *  Created on: Nov 20, 2011
 *      Author: Seungsu KIM
 */

#ifndef __GMRDYNAMICS_H__
#define __GMRDYNAMICS_H__

#include "Gaussians.h"

#define GMR_ERROR_TOLERANCE 0.001
#define INTEGRATION_L 5
#define REACHING_ITERATION_MAX 500
#define REAL_MIN (1e-30)

// GMR Dynamics
class GMRDynamics
{
private:
	Gaussians *GMM;

	double delta_t;
	double target_t;
	double current_t;

	Vector gXi;
	Vector target;
	unsigned int gDim;

public:
	GMRDynamics(int nStates, int nVar, double delta_t, const char *f_mu, const char *f_sigma, const char *f_prio );

	void initGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex);

	void   setStateTarget(Vector state, Vector target);
	void   setTarget(Vector target, double target_t= -1.0);
	Vector getTarget(void);
	double getTargetT(void);
	void   setState(Vector state);
	Vector getState(void);
	void   setCurrentTime(double current_t);
	double getCurrentTime(void);

	Vector getVelocity(Vector x);

	Vector getNextState(void);
	Vector getNextState(double lamda);
	double getReachingTime(double lamda);
};



#endif //__GMRDYNAMICS_H__
