/*
 * Copyright (c) 2020 Steven Braeger
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include "PnPDynamic.hpp"

#include<Eigen/SVD>
#include<stdexcept>
#include<Eigen/QR>
#include<iostream>


static Eigen::MatrixXd precondition_in_place(Eigen::MatrixXd& Data)
{
	size_t D=Data.rows();
	Eigen::VectorXd mu=Data.rowwise().mean();
	mu[D-1]=0.0; //mean subtraction shouldn't cancel the projective space
	Data.colwise()-=mu;
	double f=1.0/Data.topRows(D-1).norm(); //frobnius norm of the non-projective data
	Data.topRows(D-1)*=f;
	Eigen::MatrixXd P=Eigen::MatrixXd::Identity(D,D)*f; //isotropic scaling P=S T
	P.col(D-1)=-mu*f;
	P(D-1,D-1)=1.0; //enforce affine transform property
	return P;
}

Eigen::MatrixXd PnP_dynamic_in_place(Eigen::MatrixXd& Y,Eigen::MatrixXd& X,bool ignore_Yw,bool allow_last_element_zero)
{
	size_t N=Y.cols();
	if(X.cols() != N)
	{
		throw std::runtime_error("Y and X don't have the same number of columns!");
	}
	size_t L=Y.rows(),M=X.rows();
	size_t l=L-1,m=M-1;
	
	Eigen::MatrixXd Px=precondition_in_place(X);
	Eigen::MatrixXd Py=precondition_in_place(Y);

	Eigen::MatrixXd Fp;
	if(!ignore_Yw)
	{
		throw std::runtime_error("Error: TODO: Yw mode doesn't work right now");
	}
	if(ignore_Yw)  //we need to use PnP with SVD because we don't have the resulting transformed w coordinates
	{
		if(l*N < (L*M-1))
		{
			throw std::runtime_error("Not enough sample points for this dimensionality of problem");
		}
		Eigen::MatrixXd A=Eigen::MatrixXd::Zero(l*N+(allow_last_element_zero ? 0 : 1),L*M);
		Eigen::MatrixXd Ablock;
		
		for(size_t i=0;i<N;i++)
		{
			Ablock=Eigen::MatrixXd::Zero(l,L*M);
			Eigen::RowVectorXd xT=X.col(i).transpose();
			Eigen::VectorXd y=Y.col(i);
			Ablock.rightCols(M)=-y.topRows(l)*xT; //outer product for homogenous section of block
			for(size_t j=0;j<l;j++)
			{
				Ablock.block(j,j*M,1,M)=xT;  //non-homogenous section of block
			}
			A.block(i*l,0,l,L*M)=Ablock;
		}
		std::cout << A << std::endl;
		
		//std::cout << "A:\n" << A << std::endl;
		//Y=FX
		//Y'=Py Y, X'=Px X
		//F'X'=Y'
		
		//Solve for F'
		Eigen::VectorXd h;   //Ah=0
		if(allow_last_element_zero) //since the last element can be zero we have to 
		{
			//http://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD.html
			Eigen::BDCSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullV);
			h=svd.matrixV().col(A.cols()-1);//last column of V minimizes ||Ah|| s.t. ||h||=1
		}
		else
		{
			Eigen::VectorXd z=Eigen::VectorXd::Zero(A.cols());
			z[A.cols()-1]=1.0;
			A.row(A.rows()-1).array()=0.0;
			A(A.rows()-1,A.cols()-1)=1.0;
			Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
			
			h=qr.solve(z);
		}
		Fp=Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(h.data(),L,M); //decode F' from h
	}	
	else //use a least squares solution because we know the correct projections
	{
		//AX=B
		//F' X'=Y'
		//X'T F'T=Y'T
		//F'T=X'T.solve(Y'T)
		//F'=(X'T.solve(Y'T))T
		X.transposeInPlace();
		Y.transposeInPlace();
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
		Fp=qr.solve(Y);
		Fp.transposeInPlace();

	}
	//std::cout << "Fp:\n" << Fp << std::endl;
	//F'Px X=Py Y
	//inv(Py) F' Px X=Y
	//(inv(Py) F' Px) X=Y
	//F=inv(Py) F' Px
	Eigen::MatrixXd result=(Py.inverse()*Fp)*Px;
	return result;
}

