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


#include "jt_ds.h"


jt_ds::jt_ds(double dt_,int Num_C_,int Num_J_,int Num_Com_, Vector P_Joint_min_, Vector P_Joint_max_){

	dt=dt_;
	Num_C=Num_C_;
	Num_J=Num_J_;
	Num_Com=Num_Com_;
	Num_Latent=0;
	Power_of_S=2;


	Theta.Resize(Num_Com);			Theta.Zero();

	q.Resize(Num_J); 				q.Zero();
	Dq.Resize(Num_J); 				Dq.Zero();
	q_new.Resize(Num_J); 			q_new.Zero();
	q_min.Resize(Num_J); 			q_min.Zero();
	q_max.Resize(Num_J);			q_max.Zero();
	Phi.Resize(Num_Latent);			Phi.Zero();

	P_End.Resize(Num_C);			P_End.Zero();
	Target.Resize(Num_C);			Target.Zero();

	P_Matrix.Resize(Num_C,Num_C);	P_Matrix.Zero();
	Jacobian.Resize(Num_C,Num_J);	Jacobian.Zero();
	Jacobian_T.Resize(Num_J,Num_C);	Jacobian_T.Zero();

	S_q.Resize(Num_J,Num_J);		S_q.Zero();
	A_q.Resize(Num_J,Num_J);		A_q.Zero();





	if (P_Joint_min_.Size()==Num_J){
		q_min=P_Joint_min_;
	}
	else{
		cout<<"Size of P_Joint_min_ is wrong"<<endl;
	}

	if (P_Joint_max_.Size()==Num_J){
		q_max=P_Joint_max_;
	}
	else{
		cout<<"Size of P_Joint_max_ is wrong"<<endl;
	}


}

void jt_ds::initialize_P(const char  *path_){

	Matrix fMatrix; fMatrix.Resize(Num_C,Num_C);
	fMatrix.Load(path_);
	if ((fMatrix.RowSize()==Num_C)&&(fMatrix.ColumnSize()==Num_C))	{
		P_Matrix=fMatrix;
	}
	else{
		cout<<"Size of P is wrong"<<endl;
	}
	P_Matrix.Print("P_Matrix");
}
void jt_ds::initialize_A(const char  *path_){

	Matrix fMatrix; fMatrix.Resize(Num_Com*Num_J,Num_J);
	fMatrix.Load(path_);
	int j = 0;
	for(int s=0; s<Num_Com; s++ ){
		A_Matrix[s].Resize(Num_J,Num_J);
	}
	for(int s=0; s<Num_Com; s++ ){
		for(int i=0; i<Num_J; i++ ){
			A_Matrix[s].SetRow(fMatrix.GetRow(j), i);
			j++;
		}
	}
	for(int s=0; s<Num_Com; s++ ){
		A_Matrix[s].Print("A_Matrix[s]");
	}

}
void jt_ds::initialize_GMM_Latend(const char  *path_prior_,const char  *path_mu_,const char  *path_sigma_,const char  *path_W_){

	Matrix fMatrix;
	fMatrix.Load(path_W_);
	if (fMatrix.ColumnSize()==Num_J)
	{
		Num_Latent=fMatrix.RowSize();
		W_Matrix.Resize(Num_Latent,Num_J);W_Matrix.Zero();
		W_Matrix=fMatrix;
		cout<<"The size of the latent space is "<<Num_Latent<<"."<<endl;
	}
	else
	{
		cout<<"The size  of W in the latent space is not right."<<endl;
		while (ros::ok){};
	}
	W_Matrix.Print("W_Matrix");

	fMatrix.Resize(0,0);fMatrix.Zero();
	fMatrix.Load(path_prior_);
	Prior.Resize(Num_Com);			Prior.Zero();
	int nbStates=0;
	if( fMatrix(0, fMatrix.ColumnSize()-1) > 0.0 )
	{
		nbStates = fMatrix.ColumnSize();
	}
	else
	{
		nbStates = fMatrix.ColumnSize()-1;
	}

	if (Num_Com==nbStates)
	{
		for (int i=0; i<Num_Com; i++)
		{
			Prior(i)=fMatrix(0,i);
		}
	}
	else
	{
		cout<<"The number of components are wrong."<<"Num_Com= "<<Num_Com<<" nbStates="<<nbStates<<endl;
		fMatrix.Print("Prior");
		while (ros::ok){};
	}
	Prior.Print("Prior");
	for(int s=0; s<Num_Com; s++ )
	{
		Mu[s].Resize(Num_Latent);
	}
	fMatrix.Resize(0,0);fMatrix.Zero();
	fMatrix.Load(path_mu_);
	if ((fMatrix.ColumnSize()==Num_Com)&&(fMatrix.RowSize()==Num_Latent))
	{
		for(int s=0; s<Num_Com; s++ )
		{
			Mu[s].Set(fMatrix.GetColumn(s));

		}
	}
	else
	{
		cout<<"The size of matrix Mu is wrong."<<"fMatrix.ColumnSize()= "<<fMatrix.ColumnSize()<<" Num_Com= "<< Num_Com
				<<"fMatrix.RowSize()= "<<fMatrix.RowSize()<<" Num_Latent "<<Num_Latent<<endl;
		while (ros::ok){};
	}
	for(int i=0; i<Num_Com; i++ ){
		cout<<"Mu["<<i<<"]."<<endl;
		Mu[i].Print("Mu[s]");
	}

	fMatrix.Resize(0,0);fMatrix.Zero();
	fMatrix.Load(path_sigma_);

	int s=0;
	if ((fMatrix.ColumnSize()==Num_Latent)&&(fMatrix.RowSize()==Num_Com*Num_Latent))
	{
		for(int i=0; i<Num_Com; i++ ){
			Sigma[i].Resize(Num_Latent,Num_Latent);
		}
		for(int i=0; i<Num_Com; i++ ){
			for(int j=0; j<Num_Latent; j++ ){
				Sigma[i].SetRow(fMatrix.GetRow(s), j);
				s++;
			}
		}
	}
	else
	{
		cout<<"The size of matrix Sigma is wrong."<<"fMatrix.ColumnSize()= "<<fMatrix.ColumnSize()<<" Num_Latent= "<< Num_Latent
				<<"fMatrix.RowSize()= "<<fMatrix.RowSize()<<" Num_Com*Num_Latent "<<Num_Com*Num_Latent<<endl;
		while (ros::ok){};
	}
	for(int i=0; i<Num_Com; i++ ){
		cout<<"Sigma["<<i<<"]."<<endl;
		Sigma[i].Print();
	}

}

void jt_ds::Set_State(Vector P_END_,Vector P_Joints_,Matrix Jacobian_, Vector Target_){

	if (P_END_.Size()!=P_End.Size())
	{
		cout<<"Size of P_End is wrong"<<endl;
		while (ros::ok){};
	}
	P_End=P_END_;
	if (P_Joints_.Size()!=q.Size())
	{
		cout<<"Size of P_Joints_ is wrong"<<endl;
		while (ros::ok){};
	}
	for  (int i=0;i<Num_J;i++)
	{
		if ((q_min(i)<=P_Joints_(i))&&(P_Joints_(i)<=q_max(i))){
			q(i)=P_Joints_(i);
		}
		else
		{
			cout<<"We passed the joint limits. The "<<i<<" joint is outside of the defined limits. q_min(i) "<<q_min(i)
				<<" q_max(i) "<<q_max(i)<<" P_Joints_(i) "<<P_Joints_(i)<<endl;
			while (ros::ok){};
		}
	}

	if ((Jacobian_.ColumnSize()!=Jacobian.ColumnSize())&&(Jacobian_.RowSize()!=Jacobian.RowSize()))
	{
		cout<<"Size of P_Joints_ is wrong"<<endl;
		while (ros::ok){};
	}
	Jacobian=Jacobian_;
	Jacobian_T=Jacobian.Transpose();


	if (Target_.Size()==Num_C)
	{
		Target=Target_;
	}
	else
	{
		cout<<"Size of Target_ is wrong"<<endl;
		while (ros::ok){};
	}
}

void jt_ds::Update(){

	S_q=Calculate_S(q);
	A_q=Calculate_A(q);

	Dq=-S_q*A_q*S_q*Jacobian_T*P_Matrix*(P_End-Target);
/*	Jacobian_T.Print("Jacobian_T");
	P_End.Print("P_End");
	Target.Print("Target");
	P_Matrix.Print("P_Matrix");
	Dq=-A_q*Jacobian_T*P_Matrix*(P_End-Target);*/
/*	Dq.Print("Dq");*/
/*	Dq.Print("Dq");
	S_q.Print("S_q");
	A_q.Print("A_q");
	Jacobian_T.Print("Jacobian_T");
	P_Matrix.Print("P_Matrix");
	P_End.Print("P_End");
	Target.Print("Target");*/

	q_new=q+Dq*dt;
}


void jt_ds::Get_State(Vector &V_joints_,Vector &P_joints_New_){

	V_joints_=Dq;
	P_joints_New_=q_new;
}

Matrix jt_ds::Calculate_S(Vector q_){

	Matrix S_; S_.Resize(Num_J,Num_J);
	S_.Zero();

	double handle=0;
	for (int i=0;i<Num_J;i++)
	{
		handle=2*((q_(i)-q_min(i))/(q_max(i)-q_min(i)))-1;
		S_(i,i)=1-pow(handle,2*Power_of_S);
	}
//	S_.Print("S_");
	return S_;
}


Matrix jt_ds::Calculate_A(Vector q_){

	Matrix A_; A_.Resize(Num_J,Num_J);A_.Zero();

	Theta=Calculate_Theta(q_);
	for (int i=0;i<Num_Com;i++)
	{
		A_=A_+A_Matrix[i]*Theta(i);
	}
//	A_.Print("A_");
	return A_;
}


Vector jt_ds::Calculate_Theta(Vector q_)
{
	Vector Theta_;Theta_.Resize(Num_Com);Theta_.Zero();

	Phi.Zero();
	Phi=W_Matrix.Mult(q_);
//	Phi.Print("Phi");
//	W_Matrix.Print("W_Matrix");
//	Phi.Print("Phi");
	for (int i=0;i<Num_Com;i++)
	{
		Theta_(i)=Prior(i)*GaussianPDF(Phi,Mu[i],Sigma[i]);
	}
	double sum=Theta_.Sum();
	if (sum<1e-100)
	{
		for (int i=0;i<Num_Com;i++)
		{
			Theta_(i)=1.0/Num_Com;
		}
	}
	else
	{
		Theta_.Div(sum,Theta_);
	}
//	Theta_.Print("Theta_");
	return Theta_;
}


double jt_ds::GaussianPDF(Vector x_,Vector Mu_,Matrix Sigma_)
{

	double p;
	Matrix gfDiff;gfDiff.Resize(1,Num_Latent);
	Matrix gfDiff_T;gfDiff_T.Resize(Num_Latent,1);
	Matrix SigmaIIInv_;SigmaIIInv_.Resize(Num_Latent,Num_Latent);
	double detSigmaII=0;
	Matrix gfDiffp;gfDiffp.Resize(1,1);gfDiffp.Zero();

	Sigma_.Inverse(SigmaIIInv_,&detSigmaII);
	if (detSigmaII<0)
	{
		detSigmaII=0;
	}
	gfDiff.SetRow(x_ - Mu_,0);
	gfDiff_T=gfDiff.Transpose();
	gfDiffp =gfDiff*SigmaIIInv_* gfDiff_T;
	gfDiffp(0,0)=fabs(0.5*gfDiffp(0,0));

	p = exp(-gfDiffp(0,0)) / sqrt(pow(2.0*PI, Num_Latent)*( detSigmaII +1e-50));
	return p;
}


void jt_ds::Print(){
	q.Print("q");
	P_End.Print("P_End");
	Jacobian.Print("Jacobian");
	Jacobian_T.Print("Jacobian_T");
	Target.Print("Target");
}
