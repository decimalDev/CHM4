#include "Task I.h"

//из графических соображений: всего 3 корня, первый лежит в отрезке [0;0.5], второй -- на отрезке [3.5;4.5], третий -- на отрезке [4.5;5]. Положим эти точки в вектор vector<vector<double>> Starts

double f(double x)
{
	return(x + 4 * sin(x) - 1);
}

double phi(double x) //для построения метода простых итераций для первого корня воспользуемся представлением x = arcsin([1-x]/4); производная функции в правой части имеет модуль меньше 1 на большом интервале, включающем первый и второй корни уравнения, но при выбранных нами приближениях соответствующая последовательность сходится именно к первому корню
{
	return(asin((1 - x) / 4));
}

double fi(double x) //для построения метода простых итераций для второго корня воспользуемся представлением x = x + f(x) -- производная правой части имеет модуль меньше 1 в небольшой окрестности второго корня
{
	return(2 * x + 4 * sin(x) - 1);
}

double varphi(double x) //третий корень будем искать как x = 1 - 4sin(x), модуль функции в правой части меньше единицы поблизости лишь третьего из корней
{
	return(1 - 4 * sin(x));
}

double fnewton(double x)
{
	return((x + 4 * sin(x) - 1) / (1 + 4 * cos(x)));
}

void Dichotomia(vector<double>& X, double a, double b, double eps,int &iteration) //примитивная реализация метода дихотомии -- ищет корни только на тех отрезках, на концах которых функция принимает различные значения; если при делениях корень лежит в отрезке, на концах которого функция принимает одинаковое значение, то он не будет найден
{
	iteration++;
	double x = (a + b) / 2;
	if (b - a <= 2*eps || f(x) == 0)
	{
		X.push_back(x);
	}
	if(b - a > 2 * eps)
	{
		if (f(a) * f(x) <= 0)
		{
			Dichotomia(X, a, x, eps,iteration);
		}
		if (f(x) * f(b) <= 0)
		{
			Dichotomia(X, x, b, eps, iteration);
		}
	}
}

void Dichotomia_adhoc(vector<double>& X, vector<vector<double>>& I, double eps) //I -- множество отрезков (каждый представлен парой своих границ), на которых функция имеет единственный нуль
{
	int iteration = 0;
	
	for (vector<double> i : I)
	{
		double a = i[0];
		double b = i[1];
		Dichotomia(X, a, b, eps, iteration);
	}
	sort(X.begin(), X.end());
	cout << "By the dichotomia method " << X.size() << " roots were found: ";
	for (double x : X) cout << x << " ";
	cout << endl;

	cout <<" Dichotomia_adhoc number of iteration " << iteration << endl;
}

void Chordae(vector<double>& X, double a, double b, double eps, int i,int &iteration) //i = 1, если фиксируем левые концы хорд; 2, если правые
{
	iteration++;
	double c, fc, t, x, buf; //t -- предыдущая точка, x -- текущая
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
	int iterations = 0;
	do
	{
		iterations++;
		x = t - (t - c) * f(t) / (f(t) - fc);
		buf = t;
		t = x;
	} while (abs(buf - x) > eps);
	X.push_back(x);
	//cout << "Chordae number of iteration " << iterations << endl;
}

void Chordae_adhoc(vector<double>& X, vector<vector<double>>& I, double eps)
{
	int iterations = 0;
	X = {};
	cout << "To find the roots using the chordae method, choose which ends are to be fixed: press 1 for left ones, 2 for right ones." << endl;
	int Flag;
	cin >> Flag;
	for(vector<double> i: I)
	{
		double a = i[0];
		double b = i[1];
		Chordae(X, a, b, eps, Flag, iterations);
	}
	sort(X.begin(), X.end());
	cout << "By the chordae method " << X.size() << " roots were found: ";
	for (double x : X) cout << x << " ";
	cout << endl;
	cout << "Chordae number of iteration " << iterations << endl;
}

vector<double> Iterationes(vector<vector<double>>& I, double eps)
{
	int iterations = 0;
	vector<double> X;
	for (int i = 0; i < 3; i++)
	{
		double a = I[i].at(0);
		double b = I[i].at(1);
		double xold = (a + b) / 2;
		double xnew, buf;
		do
		{
			iterations++;
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
	cout << "Iterationes number of iteration " << iterations << endl;
	return X;
}

vector<double> Newton(vector<vector<double>>& I, double eps)
{
	int iterations = 0;
	vector<double> X;
	for (int i = 0; i < 3; i++)
	{
		iterations++;
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
	cout << "Newton number of iteration " << iterations << endl;
	return X;
}

vector<double> Secantes(vector<vector<double>>& I, double eps)
{
	int iterations = 0;
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
			iterations++;
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
	cout << "Secantes number of iteration " << iterations << endl;
	return X;
}