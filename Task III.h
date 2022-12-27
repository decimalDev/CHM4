#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.h"

using namespace std;

double det(Matrix& A);

vector<double> SysNewton(double(*f)(double, double), double(*g)(double,double), double (*dFdx)(double, double), double (*dFdy)(double, double), double (*dGdx)(double, double), double (*dGdy)(double, double), double x0, double y0, double eps);

double ftest(double x, double y);

double gtest(double x, double y);

double dFdx(double x, double y);

double dFdy(double x, double y);

double dGdx(double x, double y);

double dGdy(double x, double y);

double fsysnewton(double x, double y);

double gsysnewton(double x, double y);

double dFSNdx(double x, double y);

double dGSNdx(double x, double y);

double dFSNdy(double x, double y);

double dGSNdy(double x, double y);

double Numerical_dFdx(double x, double y);

double Numerical_dFdy(double x, double y);

double Numerical_dGdx(double x, double y);

double Numerical_dGdy(double x, double y);
