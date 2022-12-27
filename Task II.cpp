#include "Task II.h"

vector<double> SystemIter(double x0, double y0, double (*phi_x)(double, double), double (*phi_y)(double, double), double eps)
{
	double xold = x0;
	double yold = y0;
	double xnew, ynew, xbuf, ybuf;
	do
	{
		xnew = phi_x(xold, yold);
		ynew = phi_y(xold, yold);
		xbuf = xold;
		ybuf = yold;
		xold = xnew;
		yold = ynew;
	} while (sqrt((xnew-xbuf)*(xnew - xbuf) + (ynew-ybuf)* (ynew - ybuf)) > eps);
	cout << "The root of this system equals (" << xnew << "; " << ynew << ")^T." << endl;
	vector<double> r = { xnew, ynew };
	return(r);
}

double test_x(double x, double y)
{
	return(sqrt((x*(y+5)-1)/2));
}

double test_y(double x, double y)
{
	return(sqrt(x + 3*log(x)));
}

//ПОИСК МИНИМУМА: наша функция имеет вид F(x,y) = sin(xy) + y^2 + x^4
//необходимое условие эстремума dF/dx = 0 = dF/dy принимает вид
//ycos(xy) + 4x^3 = 0,
//xcos(x,y) + 2y = 0
//Из этой системы нужно получить некоторое УДАЧНОЕ представление вида
//x = phi_1 (x,y),
//y = phi_2 (x,y)
//Далеко не всякое такое представление нам подойдёт: чтобы итерация сходилась к искомому корню, частные производные функций phi_1 и phi_2 должны удовлетворять условиям, которых не так-то просто добиться
//Пришлось перебрать много вариантов, но в итоге сработал следующий: полагая, что x != 0 и y != 0 (это позволяет чертёж), в обоих уравнениях выразим cos(xy):
//cos(xy) = -4x^3/y
//cos(xy) = -2y/x
//Приравнивая левые части, получаем
//y^2 = 2x^4
//Предположим, что мы работаем в полуплоскости y > 0 (из графических соображений мы ищем два симметричных минимума, один из которых лежит в нужной полуплоскости), тогда
//y = sqrt(2) * x^2
//Функция phi_2 получена. Чтобы получить phi_1, подставим полученное выражение для y в первое исходное условие, сократим равенство на x^2 и получим
//x = -sqrt(2) cos(sqrt(2)*x^3)/4
//Это phi_1. Полученная пара функций даёт хорошее приближение (проверил графически). При x = 0.25 нужные достаточные условия сходимости очевидно соблюдаются (график), поэтому в качестве начальной точки положим x = 0.25, y = sqrt(2)*0.25*0.25

double F_x(double x, double y)
{
	//return(-pow(y * cos(x * y) / 4, 1 / 3));
	return(-sqrt(2)*cos(sqrt(2)*pow(x,3))/4);
}

double F_y(double x, double y)
{
	return(sqrt(2)*x*x);
}