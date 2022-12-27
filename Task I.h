#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.h"

using namespace std;

double f(double x);

double phi(double x);

double fi(double x);

double varphi(double x);

double fnewton(double x);

void Dichotomia(vector<double>& X, double a, double b, double eps);

void Dichotomia_adhoc(vector<double>& X, vector<vector<double>>& I, double eps);

void Chordae(vector<double>& X, double a, double b, double eps, int i);

void Chordae_adhoc(vector<double>& X, vector<vector<double>>& I, double eps);

vector<double> Iterationes(vector<vector<double>>& I, double eps);

vector<double> Newton(vector<vector<double>>& I, double eps);

vector<double> Secantes(vector<vector<double>>& I, double eps);


