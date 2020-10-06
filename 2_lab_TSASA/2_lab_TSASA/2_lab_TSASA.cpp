#include <iostream>
#include<cmath>
#include<vector>
#include <iomanip>

double Multimodal_Function(const double& x)
{
	//double y = pow((1 - x), 2) + exp(x);
	//double y = -sqrt(x) * sin(x) + 2;
	double y = cos(x) * tanh(x)*sin(5*x);
	return y;
}
double Unimodal_Function(const double& x)
{
	//double y = pow((1 - x), 2) + exp(x);
	//double y = -sqrt(x) * sin(x) + 2;
	double y = cos(x) * tanh(x);
	return y;
}
double Minimum_Value(std::vector<std::vector<double>>& point, double  y,double i, double j, int k)
{
	if (k == 0)
	{
		return y;
	}
	if (point[i][j]>y)
	{
		return y;
	}
	return point[i][j];
}
int main()
{
	double q = 0.005;
	std::vector<std::vector<int>> Amount_of_Points(20, std::vector<int>(10));
	std::cout << "+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+"<<std::endl;
	std::cout << "| q\\P\t|";
	for (double P = 0.9; P <1; P=P+0.01)
	{
		std::cout << " " << P << "\t|";
	}
	std::cout << std::endl;
	std::cout << "+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+" << std::endl;
	for (size_t i = 0; i <20; i=i++)
	{
		double P=0.9;
		std::cout << "| " << q << "\t|";
		for (size_t j = 0; j <10; j++)
		{
			Amount_of_Points[i][j] = ceil(log(1 - P) / log(1 - q));
			std::cout <<" "<< Amount_of_Points[i][j] << "\t|";
			P = P + 0.01;
		}
		std::cout<<std::endl;
		q = q + 0.005;
	}
	std::cout << "+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+" << std::endl;


	
	//Unimodal
	std::cout << std::endl;
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;
	std::cout << "| q\\P\t|";
	for (double P = 0.9; P < 1; P = P + 0.01)
	{
		std::cout << "     " << P << "\t|";
	}
	std::cout << std::endl;
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;
	std::vector<std::vector<double>> Unimodal_Value(20, std::vector<double>(10));
	q = 0.005;
	for (size_t i = 0; i < 20; i = i++)
	{
		std::cout << "| " << q << "\t|";
		for (size_t j = 0; j < 10; j++)
		{ 
			double y = 0;
			for (size_t k = 0; k <Amount_of_Points[i][j] ; k++)
			{
				
				y = Unimodal_Function((float)(rand() % 25001) / 10000 + 1.5);
				Unimodal_Value[i][j] = Minimum_Value(Unimodal_Value, y, i, j, k);
			}
			std::cout<< std::setprecision(4)<<"   " << Unimodal_Value[i][j] << "  \t|";
		}
		std::cout << std::endl;
		q = q + 0.005;
	}
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;


	//Multimodal
	std::cout << std::endl;
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;
	std::cout << "| q\\P\t|";
	for (double P = 0.9; P < 1; P = P + 0.01)
	{
		std::cout << "     " << P << "\t|";
	}
	std::cout << std::endl;
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;
	std::vector<std::vector<double>> Multimodal_Value(20, std::vector<double>(10));
	q = 0.005;
	for (size_t i = 0; i < 20; i = i++)
	{
		std::cout << "| " << q << "\t|";
		for (size_t j = 0; j < 10; j++)
		{
			double y = 0;
			for (size_t k = 0; k < Amount_of_Points[i][j]; k++)
			{

				y = Multimodal_Function((float)(rand() % 25001) / 10000 + 1.5);
				Multimodal_Value[i][j] = Minimum_Value(Multimodal_Value, y, i, j, k);
			}
			std::cout << std::setprecision(4) << "   " << Multimodal_Value[i][j] << "  \t|";
		}
		std::cout << std::endl;
		q = q + 0.005;
	}
	std::cout << "+-------+---------------+---------------+---------------+---------------+---------------+---------------+";
	std::cout << "--------------+----------------+---------------+---------------+" << std::endl;
	
}
