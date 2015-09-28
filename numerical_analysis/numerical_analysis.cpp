//Author: Sarah P. Flanagan
//Date: 28 September 2015
//Purpose: perform numerical analysis to check fst-heterozygosity relationship

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

using namespace std;

class population
{
public:
	int pop_size, sampled_p;
	double p, q, het;
};

int main()
{
	int end, time, i;
	int num_pops;
	double migration_rate, pbar, qbar, fst, last_fst, ht, hs_sum, last_p;
	double nc, s2, a, b, c, wc_fst, last_wc_fst;
	int r, nbar;
	bool equilibrium = false;
	default_random_engine generator;
	vector<population> pops;
	population total_pop;
	

	//initialize variables
	num_pops = r = 100;
	migration_rate = 0.01;
	pbar = 0.5;
	normal_distribution<double> norm_dist(pbar, 0.05); //arbitrary stddev
	qbar = 1 - pbar;
	last_fst = last_wc_fst = 0;
	time = 0;
	total_pop.pop_size = 10000;
	total_pop.p = 0;
	total_pop.q = 0;
	total_pop.het = 0;
	for (i = 0; i < num_pops; i++)
	{
		pops.push_back(population());
		pops[i].p = norm_dist(generator);
		pops[i].q = 1 - pops[i].p;
		pops[i].het = pops[i].p * pops[i].q * 2;
		pops[i].pop_size = total_pop.pop_size/num_pops;
		total_pop.p = total_pop.p + pops[i].p;
		total_pop.q = total_pop.q + pops[i].q;
		total_pop.het = total_pop.het + pops[i].het;
	}
	total_pop.p = total_pop.p / num_pops;
	total_pop.q = total_pop.q / num_pops;
	total_pop.het = total_pop.het / num_pops;

	while (!equilibrium)
	{
		total_pop.p = 0;
		total_pop.q = 0;
		total_pop.het = 0;
		nbar = 0;
		s2 = 0;
		for (i = 0; i < num_pops; i++)
		{
			last_p = pops[i].p;
			//Island model determines p
			//pt = pt-1(1-m)+pbarm
			pops[i].p = last_p*(1-migration_rate)+(pbar*migration_rate);
			//then drift happens
			binomial_distribution<int> distribution((2 * pops[i].pop_size), pops[i].p);
			pops[i].sampled_p = distribution(generator);
			pops[i].p = (double)pops[i].sampled_p / (2 * (double)pops[i].pop_size);
			pops[i].q = 1 - pops[i].p;
			pops[i].het = 2 * pops[i].p * pops[i].q;	
			total_pop.p = total_pop.p + pops[i].p;
			total_pop.q = total_pop.q + pops[i].q;
			total_pop.het = total_pop.het + pops[i].het;
			nbar = nbar + pops[i].pop_size;
		}
		total_pop.p = total_pop.p / num_pops;
		total_pop.q = total_pop.q / num_pops;
		total_pop.het = total_pop.het / num_pops;
		nbar = nbar / num_pops;
		//then calculate fst three ways:
		//Wright's: Fst=(Ht-Hs)/Ht aka Nei's: 1-sum(Hs)/HtN
		ht = 1-((total_pop.p*total_pop.p) + (total_pop.q*total_pop.q));
		fst = (ht - total_pop.het)/ht;
		fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
		//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
		//Fst=a/(a+b+c)
		//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
		//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
		//c=0.25hbar
		//nc=(rnbar-sum(n2/rnbar))/(r-1)
		nc = r*nbar;
		for (i = 0; i < num_pops; i++)
		{
			s2 = s2 + pops[i].pop_size*(pops[i].p - total_pop.p)*(pops[i].p - total_pop.p) / ((r - 1)*nbar);
			nc = nc - ((pops[i].pop_size*pops[i].pop_size) / (r*nbar));
		}
		nc = nc / (r - 1);
		a = total_pop.p*(1 - total_pop.p) - (((r - 1) / r)*s2) - 0.25*total_pop.het;
		a = s2 - ((1 / (nbar - 1))*a);
		a = a*(nbar / nc);
		//a = (nbar / nc)*(s2 - ((1 / (nbar - 1))*(total_pop.p*(1 - total_pop.p) - (((r - 1) / r)*s2) - 0.25*total_pop.het)));
		b = (nbar / (nbar - 1))*(total_pop.p*(1 - total_pop.p) - (((r - 1) / r)*s2) - ((2 * nbar - 1) / 4 * nbar)*total_pop.het);
		c = 0.5*total_pop.het;
		wc_fst = a / (a + b + c);

		if (last_fst == fst)
			equilibrium;
	
		last_fst = fst;
		last_wc_fst = wc_fst;
		time++;
		if (time % 100 == 0)
			cout << "\nTime: " << time << '\t' << fst << '\t' << wc_fst;
	}
	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}