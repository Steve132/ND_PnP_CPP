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

#ifndef PNP_DYNAMIC_HPP
#define PNP_DYNAMIC_HPP

#include<Eigen/Dense>

//http://www.cse.psu.edu/~rtc12/CSE486/lecture16.pdf
//Solves Y=FX where F is an arbitrary perspective matrix of size LxM and Y,X are data matrices of dimension LxN and MxN respectfully, where a (l+1)=L,(m+1)=M p
//Data is in each column where the last element of the column is the projective component (e.g. if this was OpenGL data it would be 4xN)
Eigen::MatrixXd PnP_dynamic_in_place(
	Eigen::MatrixXd& Y,	//The output data.
	Eigen::MatrixXd& X, //the input data
	bool ignore_Yw=true,//setting ignore_Yw=true uses an more accurate slower algorithm which assumes that the last row of the provided Y data is unknown (because you don't know the perspective divide results)
	bool allow_last_element_zero=false
	//if you set allow_last_element_zero to false then it uses an efficient solver based on QR decomposition but results could be wrong if last element is close to 0 
	//	(this param is only relevant if ignore_Yw=true) 
);

#endif
