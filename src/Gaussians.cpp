/*
 * Gaussians.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: Seungsu KIM
 */

#include <math.h>
#include <iostream>
#include <fstream>

#include "Gaussians.h"

using namespace std;
/*
Gaussians::Gaussians(GMMs *model)
{
	this->model.nbStates = model->nbStates;
	this->model.nbDim    = model->nbDim;

	this->model.States = (GMMState  *)malloc(model->nbStates*sizeof(GMMState) );

	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}
*/

Gaussians::Gaussians(const char *f_mu, const char *f_sigma, const char *f_prio)
{
	int s, i, j;
	int nbStates;
	int nbDim;

	Matrix fMatrix;
	fMatrix.Load(f_prio);
	if( fMatrix(0, fMatrix.ColumnSize()-1) > 0.0 )
	{
		nbStates = fMatrix.ColumnSize();
	}
	else
	{
		nbStates = fMatrix.ColumnSize()-1;
	}

	for( s=0; s<nbStates; s++ ){
		model.States[s].Prio = fMatrix(0,s);
	}

	fMatrix.Load(f_mu);
	nbDim = fMatrix.RowSize();
	model.nbStates = nbStates;
	model.nbDim    = nbDim;


	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu.Resize(nbDim);
		model.States[s].Sigma.Resize(nbDim, nbDim );
	}

	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu = fMatrix.GetColumn(s);
	}

	fMatrix.Load(f_sigma);
	j = 0;
	for( s=0; s<nbStates; s++ ){
		for( i=0; i<nbDim; i++ ){
			model.States[s].Sigma.SetRow(fMatrix.GetRow(j), i);
			j++;
		}
	}
}

Gaussians::Gaussians(int nbStates, int nbDim, const char *f_mu, const char *f_sigma, const char *f_prio)
{

	int s, i, j;

	model.nbStates = nbStates;
	model.nbDim    = nbDim;

	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu.Resize(nbDim);
		model.States[s].Sigma.Resize(nbDim, nbDim );
	}

	Matrix fMatrix(nbDim,nbStates);
	fMatrix.Load(f_mu);
	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu = fMatrix.GetColumn(s);
	}

	fMatrix.Resize(nbStates*nbDim,nbDim);
	fMatrix.Load(f_sigma);
	j = 0;
	for( s=0; s<nbStates; s++ ){
		for( i=0; i<nbDim; i++ ){
			model.States[s].Sigma.SetRow(fMatrix.GetRow(j), i);
			j++;
		}
	}

	fMatrix.Resize(1,nbStates);
	fMatrix.Load(f_prio);
	Vector fVector(nbStates);
	for( s=0; s<nbStates; s++ ){
		model.States[s].Prio = fMatrix(0,s);
	}

}

void Gaussians::setGMMs(GMMs *model)
{
	for(unsigned int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->States[s].Mu;
		this->model.States[s].Sigma = model->States[s].Sigma;
		this->model.States[s].Prio  = model->States[s].Prio;
	}
}


void Gaussians::InitFastGaussians(int first_inindex, int last_inindex)
{
	double det;
	int nbIN = last_inindex-first_inindex+1;

	for(unsigned int s=0; s<model.nbStates; s++ ){
		gmmpinv[s].MuI.Resize(nbIN);
		gmmpinv[s].SigmaII.Resize(nbIN, nbIN);
		gmmpinv[s].SigmaIIInv.Resize(nbIN, nbIN);
	}

	for(unsigned int s=0; s<model.nbStates; s++ ){
		for( int i=first_inindex; i<=last_inindex; i++ ) gmmpinv[s].MuI(i-first_inindex) = model.States[s].Mu(i);
		for( int i=first_inindex; i<=last_inindex; i++ )
			for( int j=first_inindex; j<=last_inindex; j++ )
				gmmpinv[s].SigmaII(i-first_inindex,j-first_inindex) = model.States[s].Sigma(i,j);

		gmmpinv[s].SigmaII.Inverse(gmmpinv[s].SigmaIIInv, &det);
		//gmmpinv[s].SigmaIIInv = gmmpinv[s].SigmaII.Inverse();
		//gmmpinv[s].SigmaII.Inverse(&det);
		if(det<0) det = 1e-30;
		gmmpinv[s].detSigmaII = det;

	}

	nbDimI = last_inindex - first_inindex +1;
	gfDiff.Resize(nbDimI);
	gfDiffp.Resize(nbDimI);
	gDer.Resize(nbDimI);

}

double Gaussians::GaussianPDFFast(int state, Vector x)
{
	double p;
	gfDiff  = x - gmmpinv[state].MuI;
	gfDiffp = gmmpinv[state].SigmaIIInv * gfDiff;

	p = exp(-0.5*gfDiff.Dot(gfDiffp)) / sqrt(pow(2.0*PI, nbDimI)*( gmmpinv[state].detSigmaII +1e-30));

	return p;
}

double Gaussians::GaussianProbFast(Vector x)
{
	double totalP=0;
	for(unsigned int s=0; s<model.nbStates; s++ ){
		totalP += model.States[s].Prio*GaussianPDFFast(s,x);
	}
	return totalP;
}

Vector Gaussians::GaussianDerProbFast(Vector x)
{
	gDer.Zero();
	for(unsigned int s=0; s<model.nbStates; s++ ){
		gDer += (gmmpinv[s].SigmaIIInv *(x-gmmpinv[s].MuI))*model.States[s].Prio*GaussianPDFFast(s,x);
	}
	return -gDer;
}

void Gaussians::InitFastGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex)
{
	double det;
	int nbIN  = last_inindex-first_inindex+1;
	int nbOUT = last_outindex-first_outindex+1;

	gPdf.Resize(model.nbStates);

	for(unsigned int s=0; s<model.nbStates; s++ ){
		gmmpinv[s].MuI.Resize(nbIN);
		gmmpinv[s].SigmaII.Resize(nbIN, nbIN);
		gmmpinv[s].SigmaIIInv.Resize(nbIN, nbIN);

		gmmpinv[s].muO.Resize(nbOUT);
		gmmpinv[s].SigmaIO.Resize(nbIN, nbOUT);
		gmmpinv[s].SigmaIOInv.Resize(nbOUT, nbOUT);
	}

	for(unsigned int s=0; s<model.nbStates; s++ ){
		for( int i=first_inindex; i<=last_inindex; i++ ){
			gmmpinv[s].MuI(i-first_inindex) = model.States[s].Mu(i);

			for( int j=first_inindex; j<=last_inindex; j++ ){
				gmmpinv[s].SigmaII(i-first_inindex,j-first_inindex) = model.States[s].Sigma(i,j);
			}
			for( int j=first_outindex; j<=last_outindex; j++ ){
				gmmpinv[s].SigmaIO(i-first_inindex,j-first_outindex) = model.States[s].Sigma(i,j);
			}
		}

		for( int i=first_outindex; i<=last_outindex; i++ ){
			gmmpinv[s].muO(i-first_outindex) = model.States[s].Mu(i);
		}

		gmmpinv[s].SigmaII.Inverse(gmmpinv[s].SigmaIIInv, &det);
		if(det<0) det = 1e-30;
		gmmpinv[s].detSigmaII = det;
		(gmmpinv[s].SigmaIO).Transpose().Inverse(gmmpinv[s].SigmaIOInv, &det);
	}

	nbDimI = last_inindex - first_inindex +1;
	gfDiff.Resize(nbDimI);
	gfDiffp.Resize(nbDimI);
	gDer.Resize(nbDimI);

}

void Gaussians::Regression(const Vector & indata, Vector & outdata, Matrix & derGMR)
{
	Regression(indata, outdata);
	cout << "derivative is not implemented yet " << endl;
}

void Gaussians::Regression(const Vector & indata, Vector & outdata)
{
	double pdfall;
	Vector h(model.nbStates);
	Vector r_diff(outdata.Size());

	for(unsigned int s=0; s<model.nbStates; s++){
		gPdf(s) = model.States[s].Prio*GaussianPDFFast(s, indata);
	}
	pdfall = gPdf.Sum();

	outdata.Zero();
	for(unsigned int s=0; s<model.nbStates; s++){
		//h(s) = gPdf(s)/(pdfall + 1e-30 );
		h(s) = gPdf(s)/(pdfall );
		r_diff = gmmpinv[s].SigmaIO.Transpose() * gmmpinv[s].SigmaIIInv * (indata - gmmpinv[s].MuI);

		for(unsigned int i=0; i<r_diff.Size(); i++ ){
			outdata(i) += h(s)*(r_diff(i) + gmmpinv[s].muO(i));
		}
	}
}

Vector Gaussians::Regression(const Vector & indata)
{
	Vector outdata(indata.Size());
	Regression(indata, outdata);
	return outdata;
}


/*
#include <math.h>
#include "Gaussians.h"

#include "armadillo"

using namespace arma;
using namespace std;

Gaussians::Gaussians(GMMs *model)
{
	this->model.nbStates = model->nbStates;
	this->model.nbDim    = model->nbDim;
	
	this->model.States = (GMMState  *)malloc(model->nbStates*sizeof(GMMState) );

	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}

Gaussians::Gaussians(int nbStates, int nbDim, char *f_mu, char *f_sigma, char *f_prio)
{
	
	int s, i, j;

	model.nbStates = nbStates;
	model.nbDim    = nbDim;
	model.States = (GMMState  *)malloc(nbStates*sizeof(GMMState) );

	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu       =  zeros<vec>(nbDim);
		model.States[s].Sigma    =  zeros<mat>(nbDim, nbDim );
	}

	// f_mu	
	ifstream fin;
	fin.open(f_mu);
	for( i=0; i<nbDim; i++ ){
		for( s=0; s<nbStates; s++ ){
			fin >> model.States[s].Mu(i);
		}
	}
	fin.close();
	
	// f_sigma
	fin.open(f_sigma);
	for( s=0; s<nbStates; s++ ){
		for( i=0; i<nbDim; i++ ){
			for( j=0; j<nbDim; j++ ){
				fin >> model.States[s].Sigma(i,j);
			}
		}
	}
	fin.close();

	// f_prio
	fin.open(f_prio);
	for( s=0; s<nbStates; s++ ){
		fin >>model.States[s].Prio;
	}
	fin.close();
}

void Gaussians::setGMMs(GMMs *model)
{
	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}


void Gaussians::InitFastGaussians(int first_inindex, int last_inindex)
{
	gmmpinv = (GMMStateP *)malloc(model.nbStates*sizeof(GMMStateP) );

	for(int s=0; s<model.nbStates; s++ ){
		gmmpinv[s].MuI = model.States[s].Mu.subvec(first_inindex, last_inindex);
		gmmpinv[s].SigmaII = model.States[s].Sigma.submat(first_inindex, first_inindex, last_inindex, last_inindex);
		gmmpinv[s].SigmaIIInv = pinv(gmmpinv[s].SigmaII);
		gmmpinv[s].detSigmaII = det(gmmpinv[s].SigmaII);
	}

	nbDimI = last_inindex - first_inindex +1;
	gfDiff  = zeros<vec>(nbDimI);
	gfDiffp = zeros<vec>(nbDimI);
	gDer    = zeros<vec>(nbDimI);
}

double Gaussians::GaussianPDFFast(int state, vec x)
{
	double p;
	gfDiff  = x - gmmpinv[state].MuI;
	gfDiffp = gmmpinv[state].SigmaIIInv * gfDiff;
	
	p = exp(-0.5*dot(gfDiff, gfDiffp)) / sqrt(pow(2.0*math::pi(), nbDimI)*( gmmpinv[state].detSigmaII +1e-30));

	return p;
}

double Gaussians::GaussianProbFast(vec x)
{
	double totalP=0;
	for(int s=0; s<model.nbStates; s++ ){
		totalP += model.States[s].Prio*GaussianPDFFast(s,x);
	}
	return totalP;
}

vec Gaussians::GaussianDerProbFast(vec x)
{
	gDer.zeros();
	for(int s=0; s<model.nbStates; s++ ){
		gDer += model.States[s].Prio * gmmpinv[s].SigmaIIInv *(x-gmmpinv[s].MuI)*GaussianPDFFast(s,x);
	}
	return -gDer;
}

//-------------------------------------------------------------------------------------------------------
void AllocateGMMs(GMMs *model, int nbDim, int nbStates)
{
	model->nbDim = nbDim;
	model->nbStates = nbStates;
	model->GMMState = (GMMState  *)malloc(nbStates*sizeof(GMMState) );

	for(int s=0; s<nbStates; s++ ){
		model->GMMState[s].Mu       =  zeros<vec>(nbDim);
		model->GMMState[s].Sigma    =  zeros<mat>(nbDim, nbDim );
	}
}


double GaussianPDF(vec x, vec Mu, mat Sigma)
{
	double p;
	vec diff  = x - Mu;
	vec diffp = pinv( Sigma ) * diff;
	int nbDim = x.size();

	p = exp(-0.5*dot(diff, diffp)) / sqrt(pow(2.0*math::pi(), nbDim)*( abs(det(Sigma)) +1e-30));

    if(p < 1e-30){
		return 1e-30;
    }
	else{
		return p;
	}
}

void GaussianMux(GMMs *modelK, GMMs *modelL, GMMs *modelOut)
{
	int k,l,j;
	int K = modelK->nbStates;
	int L = modelL->nbStates;
	int J = K*L;

	//modelOut->nbDim = modelK->nbDim;
	//modelOut->nbStates = J;
	//modelOut->GMMState = (GMMState *)malloc(J*sizeof(GMMState) );

	j=0;
	for(k=0; k<K; k++){
		for(l=0; l<L; l++){
			modelOut->GMMState[j].Sigma = pinv( pinv(modelK->GMMState[k].Sigma) + pinv( modelL->GMMState[l].Sigma) );
			modelOut->GMMState[j].Mu    = modelOut->GMMState[j].Sigma *( pinv(modelK->GMMState[k].Sigma) * modelK->GMMState[k].Mu + pinv(modelL->GMMState[l].Sigma) * modelL->GMMState[l].Mu );
			modelOut->GMMState[j].Prio  = modelK->GMMState[k].Prio* modelL->GMMState[l].Prio * GaussianPDF( modelK->GMMState[k].Mu, modelL->GMMState[l].Mu, modelK->GMMState[k].Sigma + modelL->GMMState[l].Sigma);
			j++;
		}		
	}
}

void GaussianRotate(GMMs *model, arma::vec P, arma::mat R, GMMs *modelOut)
{
	for(int s=0; s<model->nbStates; s++){
		modelOut->GMMState[s].Mu    = R*model->GMMState[s].Mu + P;
		modelOut->GMMState[s].Sigma = R*model->GMMState[s].Sigma*trans(R);
	}
	//modelOut->nbDim = model->nbDim;
	//modelOut->nbStates = model->nbStates;
}
*/
