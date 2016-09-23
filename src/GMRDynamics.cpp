/*
 * GMRDynamics.cpp
 *
 *  Created on: Nov 20, 2011
 *      Author: Seungsu KIM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GMRDynamics.h"


GMRDynamics::GMRDynamics(int nStates, int nVar, double delta_t, const char *f_mu, const char *f_sigma, const char *f_prio )
{
	this->delta_t = delta_t;
	GMM = new Gaussians(nStates, nVar, f_mu, f_sigma, f_prio);
}

void GMRDynamics::initGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex)
{
	GMM->InitFastGMR(first_inindex, last_inindex, first_outindex, last_outindex);

	gDim = last_inindex- first_inindex+1;
	if( gDim != ( last_outindex - first_outindex+1 ) ){
		cout << "dynamics dimension is not matching" << endl;
	}

	gXi.Resize(gDim);
	target.Resize(gDim);

	gXi.Zero();
	target.Zero();
}

void GMRDynamics::setStateTarget(Vector state, Vector target)
{
	setTarget(target);
	setState(state);
}

void GMRDynamics::setTarget(Vector target, double target_t)
{
	this->target_t = target_t;

	//gXi += (this->target - target);
	this->target = target;
}

Vector GMRDynamics::getTarget(void)
{
	return target;
}

double GMRDynamics::getTargetT(void)
{
	return target_t;
}

void GMRDynamics::setState(Vector state)
{
	gXi = state;
}

Vector GMRDynamics::getState(void)
{
	return gXi;
}

void GMRDynamics::setCurrentTime(double current_t)
{
    this->current_t = current_t;
}

double GMRDynamics::getCurrentTime(void)
{
    return current_t;
}

Vector GMRDynamics::getVelocity(Vector x)
{
	return GMM->Regression(x);
}


Vector GMRDynamics::getNextState(void)
{
	return getNextState(1.0);
}

Vector GMRDynamics::getNextState(double lamda)
{
	// target time
	target_t -= (delta_t*lamda);

	gXi += (getVelocity(gXi-target) *(delta_t*lamda));

	return gXi;
}

double GMRDynamics::getReachingTime(double lamda)
{
	unsigned int frame=0;
	unsigned int li=0;
	Vector xi(3);
	xi.Set(gXi);

	for(frame=0; frame<REACHING_ITERATION_MAX; frame++){
		for(li=0; li<INTEGRATION_L; li++){
			xi += getVelocity(xi-target)*delta_t/(double)INTEGRATION_L*lamda;

			if( (xi-target).Norm() < GMR_ERROR_TOLERANCE ){
				return (double)(frame*INTEGRATION_L +li)*delta_t/(double)INTEGRATION_L;
			}
		}
	}
	return (double)(frame*INTEGRATION_L +li)*delta_t/(double)INTEGRATION_L;
}

