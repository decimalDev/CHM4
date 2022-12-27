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

//����� ��������: ���� ������� ����� ��� F(x,y) = sin(xy) + y^2 + x^4
//����������� ������� ��������� dF/dx = 0 = dF/dy ��������� ���
//ycos(xy) + 4x^3 = 0,
//xcos(x,y) + 2y = 0
//�� ���� ������� ����� �������� ��������� ������� ������������� ����
//x = phi_1 (x,y),
//y = phi_2 (x,y)
//������ �� ������ ����� ������������� ��� �������: ����� �������� ��������� � �������� �����, ������� ����������� ������� phi_1 � phi_2 ������ ������������� ��������, ������� �� ���-�� ������ ��������
//�������� ��������� ����� ���������, �� � ����� �������� ���������: �������, ��� x != 0 � y != 0 (��� ��������� �����), � ����� ���������� ������� cos(xy):
//cos(xy) = -4x^3/y
//cos(xy) = -2y/x
//����������� ����� �����, ��������
//y^2 = 2x^4
//�����������, ��� �� �������� � ������������� y > 0 (�� ����������� ����������� �� ���� ��� ������������ ��������, ���� �� ������� ����� � ������ �������������), �����
//y = sqrt(2) * x^2
//������� phi_2 ��������. ����� �������� phi_1, ��������� ���������� ��������� ��� y � ������ �������� �������, �������� ��������� �� x^2 � �������
//x = -sqrt(2) cos(sqrt(2)*x^3)/4
//��� phi_1. ���������� ���� ������� ��� ������� ����������� (�������� ����������). ��� x = 0.25 ������ ����������� ������� ���������� �������� ����������� (������), ������� � �������� ��������� ����� ������� x = 0.25, y = sqrt(2)*0.25*0.25

double F_x(double x, double y)
{
	//return(-pow(y * cos(x * y) / 4, 1 / 3));
	return(-sqrt(2)*cos(sqrt(2)*pow(x,3))/4);
}

double F_y(double x, double y)
{
	return(sqrt(2)*x*x);
}