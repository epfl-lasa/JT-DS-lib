
/*
 * Copyright (C) 2016 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author: Sina Mirrazavi
 * email:   sina.mirrazavi@epfl.ch
 * website: lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#include "MathLib/MathLib.h"
#include "Gaussians.h"
#include "ros/ros.h"
#include <math.h>

#define GAUSSIAN_MAXIMUM_NUMBER 50

Matrix 			A_Matrix[GAUSSIAN_MAXIMUM_NUMBER];
Vector			Prior;
Vector			Mu[GAUSSIAN_MAXIMUM_NUMBER];
Matrix			Sigma[GAUSSIAN_MAXIMUM_NUMBER];


using namespace MathLib;

class jt_ds
{
public:

	jt_ds(double dt_,int Num_C_,int Num_J_,int Num_Com_, Vector P_Joint_min_, Vector P_Joint_max_);
	void 			initialize_A(const char  *path_);
	void 			initialize_P(const char  *path_);
	void 			initialize_GMM_Latend(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_,const char  *path_W_);
	void 			Set_State(Vector P_END_,Vector P_Joints_,Matrix Jacobian_, Vector Target_);
	void 			Update();
	void 			Get_State(Vector &V_joints_,Vector &P_joints_New_);
	void 			Print();

private:

	Matrix 			Calculate_S(Vector q_);
	Matrix			Calculate_A(Vector q_);
	Vector			Calculate_Theta(Vector q_);
	double 			GaussianPDF( Vector x_,Vector Mu_,Matrix Sigma_);


	int				Num_J;
	int				Num_Latent;
	int				Num_C;
	int				Num_Com;
	int				Power_of_S;

	double 			dt;

	Vector 			q;
	Vector 			Dq;
	Vector			q_new;
	Vector 			q_min;
	Vector 			q_max;
	Vector 			P_End;
	Vector			Theta;

	Vector			Target;

	Matrix 			P_Matrix;
	Matrix 			W_Matrix;
	Matrix 			Jacobian;
	Matrix 			Jacobian_T;
	Matrix 			S_q;
	Matrix			A_q;

	Vector			Phi;


};
