#include "Task III.h"
#define DELTA 0.000001

double det(Matrix& A) //чтобы считать определитель
{
	return(A.getelement(1, 1) * A.getelement(2, 2) - A.getelement(1, 2) * A.getelement(2, 1));
}

vector<double> SysNewton(double(*f)(double, double), double(*g)(double, double), double(*dFdx)(double, double), double(*dFdy)(double, double), double(*dGdx)(double, double), double(*dGdy)(double, double), double x0, double y0, double eps) //метод Ньютона для точных производных (вектор Delta находим методом Крамера)
{
	vector<vector<double>> W = { {dFdx(x0, y0), dFdy(x0,y0)}, {dGdx(x0, y0), dGdy(x0,y0)} };//по крамеру и по Ньютону-Рафсону
	Matrix M(W);
	vector<double> y = { f(x0, y0), g(x0, y0) };
	Vector Y(y);
	double xold = x0;
	double yold = y0;
	Vector Delta(2);
	double xnew, xbuf, ynew, ybuf;
	do
	{
		double d = det(M);
		//cout << "det = " << d << endl;
		//cout << M.getelement(1, 1) <<"*" << M.getelement(2, 2) << " - " << M.getelement(1, 2) <<"*"<< M.getelement(2, 1) << endl;
		double xi = -Y.getelement(1) * M.getelement(2, 2) + Y.getelement(2) * M.getelement(1, 2);
		double eta = -Y.getelement(2) * M.getelement(1, 1) + Y.getelement(2) * M.getelement(2, 1);
		Delta.setelement(1, xi / d);
		Delta.setelement(2, eta / d);
		xnew = xold + Delta.getelement(1);
		ynew = yold + Delta.getelement(2);
		xbuf = xold;
		ybuf = yold;
		xold = xnew;
		yold = ynew;
		M.setelement(1, 1, dFdx(xnew, ynew));
		M.setelement(1, 2, dFdy(xnew, ynew));
		M.setelement(2, 1, dGdx(xnew, ynew));
		M.setelement(2, 2, dGdy(xnew, ynew));
		Y.setelement(1, f(xnew, ynew));
		Y.setelement(2, g(xnew, ynew));
		cout << xnew << " and " << ynew << endl;
	} while (sqrt((xnew - xbuf) * (xnew - xbuf) + (ynew - ybuf) * (ynew - ybuf)) > eps);
	cout << "The root found by the Newton method: (" << xnew << "; " << ynew << ")^T." << endl;
	vector<double> R = { xnew, ynew };

	return R;
}

double ftest(double x, double y)
{
	return(sin(x - 0.6) - y - 1.6);
}

double gtest(double x, double y)
{
	return(3 * x - cos(y) - 0.9);
}

double dFdx(double x, double y)
{
	return(cos(x - 0.6));
}

double dFdy(double x, double y)
{
	return(-1);
}

double dGdx(double x, double y)
{
	return(3);
}

double dGdy(double x, double y)
{
	return(sin(y));
}

double fsysnewton(double x, double y)
{
	//return(x + sqrt(2) * cos(sqrt(2) * pow(x, 3) / 4));
	return(y * cos(x * y) + 4 * x * x * x);
}

double gsysnewton(double x, double y)
{
	//return(y - sqrt(2) * x * x);
	return(x * cos(x * y) + 2 * y);
}

double dFSNdx(double x, double y)
{
	//return(1 - (3.0 / 2.0)*x*x*sin(sqrt(2)*x*x*x/4));
	return(12 * x * x - y * y * sin(x * y));
}

double dGSNdx(double x, double y)
{
	//return(-2 * sqrt(2) * x);
	return(- x * y * sin(x * y));
}

double dFSNdy(double x, double y)
{
	//return 0;
	return (cos(x * y) - y * sin(x * y));
}

double dGSNdy(double x, double y)
{
	//return 1;
	return(2 - x * x * sin(x * y));
}

double Numerical_dFdx(double x,double y) {
	double dx = 1;
	//double dy = 1;
	double d1fx = fsysnewton(x + dx, y) - fsysnewton(x - dx, y); d1fx /= (2 * dx);
	dx /= 2;
	double d2fx = fsysnewton(x + dx, y) - fsysnewton(x - dx, y); d2fx /= (2 * dx);
	d1fx = abs(d1fx - d2fx);
	while (d1fx>DELTA) {
		d1fx = d2fx;
		dx /= 2;
		d2fx = fsysnewton(x + dx, y) - fsysnewton(x - dx, y); d2fx /= (2 * dx);
		d1fx = abs(d1fx - d2fx);
		//cout << "d2fx = " << d2fx << endl;
	}
	return d2fx;
}

double Numerical_dFdy(double x, double y) {
	//double dx = 1;
	double dy = 1;
	double d1fy = fsysnewton(x , y + dy) - fsysnewton(x, y - dy); d1fy /= (2 * dy);
	dy /= 2;
	double d2fy = fsysnewton(x, y + dy) - fsysnewton(x, y - dy); d2fy /= (2 * dy);
	d1fy = abs(d1fy - d2fy);
	while (d1fy > DELTA) {
		d1fy = d2fy;
		dy /= 2;
		d2fy = fsysnewton(x, y + dy) - fsysnewton(x, y - dy); d2fy /= (2 * dy);
		d1fy = abs(d1fy - d2fy);
		//cout << "d2fy = " << d2fy << endl;
	}
	return d2fy;
}

double Numerical_dGdx(double x, double y) {
	double dx = 1;
	//double dy = 1;
	double d1fx = gsysnewton(x + dx, y) - gsysnewton(x - dx, y); d1fx /= (2 * dx);
	dx /= 2;
	double d2fx = gsysnewton(x + dx, y) - gsysnewton(x - dx, y); d2fx /= (2 * dx);
	d1fx = abs(d1fx - d2fx);
	while (d1fx > DELTA) {
		d1fx = d2fx;
		dx /= 2;
		d2fx = gsysnewton(x + dx, y) - gsysnewton(x - dx, y); d2fx /= (2 * dx);
		d1fx = abs(d1fx - d2fx);
		//cout << "d2gx = " << d2fx << endl;
	}
	return d2fx;
}

double Numerical_dGdy(double x, double y) {
	//double dx = 1;
	double dy = 1;
	double d1fy = gsysnewton(x , y + dy) - gsysnewton(x, y - dy); d1fy /= (2 * dy);
	dy /= 2;
	double d2fy = gsysnewton(x, y + dy) - gsysnewton(x, y - dy); d2fy /= (2 * dy);
	d1fy = abs(d1fy - d2fy);
	while (d1fy > DELTA) {
		d1fy = d2fy;
		dy /= 2;
		d2fy = gsysnewton(x, y + dy) - gsysnewton(x, y - dy); d2fy /= (2 * dy);
		d1fy = abs(d1fy - d2fy);
		//cout << "d2gy = " << d2fy << endl;
	}
	return d2fy;
}

