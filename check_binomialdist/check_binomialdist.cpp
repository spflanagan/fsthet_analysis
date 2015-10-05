//Author: Sarah P. Flanagan
//Date: 29 September 2015
//Purpose: check the <random> binomial_distribution program
//binomial distribution: Pr(k|n,p)=(n choose k)p^k(1-p)^(n-k)
//mean = np, variance = np(1-p)

#include <random>
#include <iostream>
#include <vector>

using namespace std;

int main()
{
	int end, num, n;
	double p, mean, var;
	vector<int> Prk;
	random_device rd;
	default_random_engine generator(rd());

	//set variables
	num = 10000;
	n = 100;
	p = 0.05;
	cout << "n: " << n << "\tp: " << p << "\tNumReps: " << num << '\n';

	mean = var = 0;
	for (int i = 0; i < num; i++)
	{
		binomial_distribution<int> distribution(n, p);
		Prk.push_back(distribution(generator));
		mean = mean + Prk[i];
	}
	mean = mean / num;

	for (int i = 0; i < num; i++)
		var = var + (Prk[i] - mean)*(Prk[i] - mean);
	var = var / num;

	cout << "\nExpected mean: " << n*p << "\tExpected Variance: " << n*p*(1 - p);
	cout << "\nActual mean: " << mean << "\tActual Variance: " << var;
	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}