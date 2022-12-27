#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.h"

using namespace std;

vector<double> SystemIter(double x0, double y0, double (*phi_x)(double, double), double (*phi_y)(double, double), double eps);

double test_x(double x, double y);

double test_y(double x, double y);

double F_x(double x, double y);

double F_y(double x, double y);
