/*
 * Gaussians.h
 *
 *  Created on: Nov 19, 2011
 *      Author: Seungsu KIM
 */

#ifndef __GAUSSIANSM_H__
#define __GAUSSIANSM_H__

#include "MathLib.h"

#define GAUSSIAN_MAXIMUM_NUMBER 50

using namespace MathLib;

struct GMMState {
	Vector Mu;
	Matrix Sigma;
	double Prio;
};

struct GMMStateP {
	Vector MuI;
	Matrix SigmaII;
	Matrix SigmaIIInv;
	double detSigmaII;

	// for GMR
	Vector muO;
	Matrix SigmaIO;
	Matrix SigmaIOInv;
};

struct GMMs{
	unsigned int nbStates;
	unsigned int nbDim;

	GMMState  States[GAUSSIAN_MAXIMUM_NUMBER];
};

class Gaussians
{
private: 
	GMMStateP gmmpinv[GAUSSIAN_MAXIMUM_NUMBER];

public:
	GMMs model;

	Gaussians(const char *f_mu, const char *f_sigma, const char *f_prio);
	Gaussians(int nbStates, int nbDim, const char *f_mu, const char *f_sigma, const char *f_prio);
	Gaussians(GMMs *model);	

	void setGMMs(GMMs *model);

	// For fast computation of GaussianPDF
	Vector gfDiff, gfDiffp;
	Vector gDer;
	Vector gPdf;
	int nbDimI;


	void InitFastGaussians(int first_inindex, int last_inindex);
	double GaussianPDFFast(int state, Vector x);
	double GaussianProbFast(Vector x);
	Vector GaussianDerProbFast(Vector x);

	void InitFastGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex);
	void Regression(const Vector & indata, Vector & outdata, Matrix & derGMR);
	void Regression(const Vector & indata, Vector & outdata);
	Vector Regression(const Vector & indata);

};
/*
void GaussianMux(GMMs *modelK, GMMs *modelL, GMMs *modelOut);
void GaussianRotate(GMMs *model, Vector P, Matrix R, GMMs *modelOut);
*/

#endif //__GAUSSIANS_H__
