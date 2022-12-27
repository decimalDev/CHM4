#include "Task I.h"

//�� ����������� �����������: ����� 3 �����, ������ ����� � ������� [0;0.5], ������ -- �� ������� [3.5;4.5], ������ -- �� ������� [4.5;5]. ������� ��� ����� � ������ vector<vector<double>> Starts

double f(double x)
{
	return(x + 4 * sin(x) - 1);
}

double phi(double x) //��� ���������� ������ ������� �������� ��� ������� ����� ������������� �������������� x = arcsin([1-x]/4); ����������� ������� � ������ ����� ����� ������ ������ 1 �� ������� ���������, ���������� ������ � ������ ����� ���������, �� ��� ��������� ���� ������������ ��������������� ������������������ �������� ������ � ������� �����
{
	return(asin((1 - x) / 4));
}

double fi(double x) //��� ���������� ������ ������� �������� ��� ������� ����� ������������� �������������� x = x + f(x) -- ����������� ������ ����� ����� ������ ������ 1 � ��������� ����������� ������� �����
{
	return(2 * x + 4 * sin(x) - 1);
}

double varphi(double x) //������ ������ ����� ������ ��� x = 1 - 4sin(x), ������ ������� � ������ ����� ������ ������� ���������� ���� �������� �� ������
{
	return(1 - 4 * sin(x));
}

double fnewton(double x)
{
	return((x + 4 * sin(x) - 1) / (1 + 4 * cos(x)));
}

void Dichotomia(vector<double>& X, double a, double b, double eps) //����������� ���������� ������ ��������� -- ���� ����� ������ �� ��� ��������, �� ������ ������� ������� ��������� ��������� ��������; ���� ��� �������� ������ ����� � �������, �� ������ �������� ������� ��������� ���������� ��������, �� �� �� ����� ������
{
	double x = (a + b) / 2;
	if (b - a <= 2*eps || f(x) == 0)
	{
		X.push_back(x);
	}
	if(b - a > 2 * eps)
	{
		if (f(a) * f(x) <= 0)
		{
			Dichotomia(X, a, x, eps);
		}
		if (f(x) * f(b) <= 0)
		{
			Dichotomia(X, x, b, eps);
		}
	}
}

void Dichotomia_adhoc(vector<double>& X, vector<vector<double>>& I, double eps) //I -- ��������� �������� (������ ����������� ����� ����� ������), �� ������� ������� ����� ������������ ����
{
	for (vector<double> i : I)
	{
		double a = i[0];
		double b = i[1];
		Dichotomia(X, a, b, eps);
	}
	sort(X.begin(), X.end());
	cout << "By the dichotomia method " << X.size() << " roots were found: ";
	for (double x : X) cout << x << " ";
	cout << endl;
}

void Chordae(vector<double>& X, double a, double b, double eps, int i) //i = 1, ���� ��������� ����� ����� ����; 2, ���� ������
{
	double c, fc, t, x, buf; //t -- ���������� �����, x -- �������
	if (i == 1)
	{
		t = a;
		c = b;
		fc = f(b);
	}
	else
	{
		t = b;
		c = a;
		fc = f(a);
	}
	do
	{
		x = t - (t - c) * f(t) / (f(t) - fc);
		buf = t;
		t = x;
	} while (abs(buf - x) > eps);
	X.push_back(x);
}

void Chordae_adhoc(vector<double>& X, vector<vector<double>>& I, double eps)
{
	X = {};
	cout << "To find the roots using the chordae method, choose which ends are to be fixed: press 1 for left ones, 2 for right ones." << endl;
	int Flag;
	cin >> Flag;
	for(vector<double> i: I)
	{
		double a = i[0];
		double b = i[1];
		Chordae(X, a, b, eps, Flag);
	}
	sort(X.begin(), X.end());
	cout << "By the chordae method " << X.size() << " roots were found: ";
	for (double x : X) cout << x << " ";
	cout << endl;
}

vector<double> Iterationes(vector<vector<double>>& I, double eps)
{
	vector<double> X;
	for (int i = 0; i < 3; i++)
	{
		double a = I[i].at(0);
		double b = I[i].at(1);
		double xold = (a + b) / 2;
		double xnew, buf;
		do
		{
			if (i == 0) xnew = phi(xold);
			if (i == 1) xnew = fi(xold);
			if (i == 2) xnew = varphi(xold);
			buf = xold;
			xold = xnew;
		} while (abs(xnew - buf) > eps);
		X.push_back(xnew);
	}
	cout << "Roots found by the simple iterations method: ";
	for (double x : X) cout << x << " ";
	cout << endl;
	return X;
}

vector<double> Newton(vector<vector<double>>& I, double eps)
{
	vector<double> X;
	for (int i = 0; i < 3; i++)
	{
		double a = I[i].at(0);
		double b = I[i].at(1);
		double xold = (a + b) / 2;
		double xnew, buf;
		do
		{
			xnew = xold - fnewton(xold);
			buf = xold;
			xold = xnew;
		} while (abs(xnew - buf) > eps);
		X.push_back(xnew);
	}
	cout << "Roots found by the Newton method: ";
	for (double x : X) cout << x << " ";
	cout << endl;
	return X;
}

vector<double> Secantes(vector<vector<double>>& I, double eps)
{
	vector<double> X;
	for (int i = 0; i < 3; i++)
	{
		double a = I[i].at(0);
		double b = I[i].at(1);
		double x_penult = b;
		double x_ult = (a + b) / 2;
		double xnew, buf;
		do
		{
			xnew = x_ult - (x_ult - x_penult)*f(x_ult)/(f(x_ult)-f(x_penult));
			x_penult = x_ult;
			buf = x_ult;
			x_ult = xnew;
		} while (abs(xnew - buf) > eps);
		X.push_back(xnew);
	}
	cout << "Roots found by the secants method: ";
	for (double x : X) cout << x << " ";
	cout << endl;
	return X;
}