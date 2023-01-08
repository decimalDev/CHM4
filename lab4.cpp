#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.h"
#include "Task I.h"
#include "Task II.h"
#include "Task III.h"

using namespace std;

int main()
{
    
    cout.precision(7);
    cout.fixed;
    
    vector<vector<double>> I = { {0, 0.5}, {3.85, 4.0}, {4.5, 4.95} };
    double pi = 2 * acos(0.0);
    double eps = pow(10.0, -6.0);

    //TASK I
    vector<double> X;
    Dichotomia_adhoc(X, I, eps);
    Chordae_adhoc(X, I, eps);
    Iterationes(I, eps);
    Newton(I, eps);
    Secantes(I, eps);
    /*
    //TASK II
    SystemIter(3.5, 2.2, test_x, test_y, 0.001); //вроде работает!
    vector<double> R = SystemIter(0.25, sqrt(2) * 0.25 * 0.25, F_x, F_y, pow(10, -5));
    //тут я проверяю, действительно ли мы близко подошли к максимуму
    //double x = R[0];
    //double y = R[1];
    //cout << sin(x * y) + y * y + pow(x, 4);
    */
    //TASK III
    //SysNewton(ftest, gtest, dFdx, dFdy, dGdx, dGdy, 0.0, -2.0, 0.001); //ВРОДЕ РАБОТАЕТ
    //удивительное дело: метод Ньютона позволяет решить пункт (2.2) безо всяких ухищрений, напрямую из необходимого условия экстремума!
    SysNewton(fsysnewton, gsysnewton, dFSNdx, dFSNdy, dGSNdx, dGSNdy, -0.25, 0.25, pow(10, -5));

    SysNewton(fsysnewton, gsysnewton, Numerical_dFdx, Numerical_dFdy, Numerical_dGdx, Numerical_dGdy, -0.25, 0.25, pow(10, -5));
    return 0;
}
