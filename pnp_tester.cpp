#include<iostream>
#include "PnPDynamic.hpp"
using namespace std;

int main(int argc,char** argv)
{
	Eigen::MatrixXd X(3,4),Y(3,4);
	
	X.col(0)=Eigen::Vector3d(1.0,1.0,1.0);
	X.col(1)=Eigen::Vector3d(2.0,2.0,1.0);
	X.col(2)=Eigen::Vector3d(1.0,2.0,1.0);
	X.col(3)=Eigen::Vector3d(3.0,4.0,1.0);
	
	Eigen::MatrixXd Fgt=Eigen::MatrixXd::Random(3,3);//Establish a ground truth projection matrix
	Fgt/=Fgt(2,2);
	
	Y=Fgt*X;
	Eigen::MatrixXd Yno_w=Y.array().rowwise()/Y.row(2).array();

	std::cout << "X:\n" << X <<std::endl;
	std::cout << "Y:" << Y <<std::endl;
	std::cout << "Yno_w:\n" << Yno_w <<std::endl;
	std::cout << "Fgt:\n" << Fgt << std::endl;
	
	
	Eigen::MatrixXd Xc,Yc;
	/*
	Xc=X;Yc=Yno_w;
	Eigen::MatrixXd Fqr=PnP_dynamic_in_place(Yc,Xc,true,false);
	Fqr/=Fqr(2,2); //normalize for scalar invariance for check
	std::cout << "Fqr error:\n" << (Fqr-Fgt).norm() << std::endl;*/
	
	Xc=X;Yc=Yno_w;
	Eigen::MatrixXd Fsvd=PnP_dynamic_in_place(Yc,Xc,true,true);
	Fsvd/=Fsvd(2,2); //normalize for scalar invariance for accuracy check
	std::cout << "Fsvd error:\n" << (Fsvd-Fgt).norm() << std::endl;
	
	Xc=X;Yc=Y;
	Eigen::MatrixXd Flsq=PnP_dynamic_in_place(Yc,Xc,false); //THIS MODE DOESN"T WORK
	Flsq/=Flsq(2,2); //normalize for scalar invariance for accuracy check
	std::cout << "Flsq error:\n" << (Flsq-Fgt).norm() << std::endl;
	return 0;
}
