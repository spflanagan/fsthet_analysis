//Author: Sarah P. Flanagan
//Date: 28 September 2015
//Purpose: perform numerical analysis to check fst-heterozygosity relationship

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <sstream>

using namespace std;

class population
{
public:
	int pop_size, sampled_p;
	double p, q, het;
};

int main(int argc, char* argv[])
{
	int end, time, i;
	int t, tt;
	int num_pops, num_reps, num_per_p, num_gens, num_sampled;
	double migration_rate, pbar, qbar, fst, ht, last_p;
	double nc, s2, a, b, c, wc_fst, pq, rs, wc_fst8, varp;
	int r, nbar, pbar_index;
	bool interactivemode = false;
	random_device rd;
	default_random_engine generator(rd());
	vector<population> pops;
	population total_pop;
	ofstream allele_freqs, output, sample_out;
	string base_name, query, tempstring1, tempstring2, output_name, sample_out_name;

	base_name = "num.analysis.";
	total_pop.pop_size = 10000;
	migration_rate = 0.01;
	num_reps = 10000;
	num_pops = 10;
	num_sampled = num_pops;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nnumerical_analysis:\n";
			cout << "Runs numerical analysis to generate Fst-Het patterns\n";
			cout << "-o:\tbase file name (include path). Example: N1000_s10_\n";
			cout << "-n:\toverall population size\n";
			cout << "-p:\tnumber of subpopulations/demes\n";
			cout << "-r:\tnumber of reps\n";
			cout << "-m:\tmigration rate\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			cout << "Input integer to quit.\n";
			cin >> end;
			return 0;
		}
		interactivemode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nnumerical_analysis:\n";
			cout << "Runs numerical analysis to generate Fst-Het patterns\n";
			cout << "-o:\tbase file name (include path). Example: N1000_s10_\n";
			cout << "-n:\toverall population size\n";
			cout << "-s:\tnumber of subpopulations/demes\n";
			cout << "-r:\tnumber of reps\n";
			cout << "-m:\tmigration rate\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-o")
			base_name = tempstring2;
		if (tempstring1 == "-n")
			total_pop.pop_size = atoi(tempstring2.c_str());
		if (tempstring1 == "-m")
			migration_rate = atof(tempstring2.c_str());
		if (tempstring1 == "-s")
			num_pops = atoi(tempstring2.c_str());
		if (tempstring1 == "-r")
			num_reps = atoi(tempstring2.c_str());
		if (tempstring1 == "-s")
			num_sampled = atoi(tempstring2.c_str());
	}

	if (interactivemode)
	{
		cout << "\nProvide base file name (include path). Example: N1000_s10_\n";
		cin >> base_name;
		cout << "\nProvide overall population size\n";
		cin >> total_pop.pop_size;
		cout << "\nProvide the number of subpopulations/demes\n";
		cin >> num_pops;
		cout << "\nProvide the number of reps\n";
		cin >> num_reps;
		cout << "\nProvide migration rate\n";
		cin >> migration_rate;
		cout << "\nProvide the number of populations to sample\n";
	}
	else
	{
		if (base_name == "num.analysis.")
			cout << "\nWARNING: base file name is set to default: num.analysis.";
		if (total_pop.pop_size == 10000)
			cout << "\nWARNING: overall population size is set to default: 10000";
		if (migration_rate == 0.01)
			cout << "\nWARNING: migration rate is set to default: 0.01";
		if (num_reps == 10000)
			cout << "\nWARNING: number of reps is set to default: 10000";
		if (num_pops == 10)
			cout << "\nWARNING: number of subpopulations is set to default: 10";
	}

	//initialize variables
	for (i = 0; i < num_pops; i++)
	{
		pops.push_back(population());
	}
	sample_out_name = base_name + "sampledpops.txt";
	sample_out.open(sample_out_name);
	sample_out << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	output_name = base_name + "output.txt";
	output.open(output_name);
	output << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	stringstream allele_freqs_name;
	allele_freqs_name << base_name << "freqs.txt";
	allele_freqs.open(allele_freqs_name.str());
	r = num_pops;
	num_gens = 5000;
	vector<population> sampled_pops;
	int num_sampled_inds = 50;
	double allele1, allele2, total_p, total_het;
	int n_samp_per_pop = (double)num_sampled / (double)num_pops;
	int samp_index = 0;
	for (t = 0; t < num_reps; t++)
	{
		//set migtrant p-value
		uniform_real_distribution <double> unidist(0.05, 0.95);
		pbar = unidist(generator);
		qbar = 1 - pbar;
		
		for (tt = 0; tt < num_gens; tt++)
		{
			total_pop.p = 0;
			total_pop.q = 0;
			total_pop.het = 0;
			if (tt == 0)
			{
				for (i = 0; i < num_pops; i++)
				{
					pops[i].p = 0.5;
					pops[i].q = 1 - pops[i].p;
					pops[i].het = pops[i].p * pops[i].q * 2;
					pops[i].pop_size = total_pop.pop_size / num_pops;
					total_pop.p = total_pop.p + pops[i].p;
					total_pop.q = total_pop.q + pops[i].q;
					total_pop.het = total_pop.het + pops[i].het;
				}
				total_pop.p = total_pop.p / num_pops;
				total_pop.q = total_pop.q / num_pops;
				total_pop.het = total_pop.het / num_pops;
			}
			nbar = 0;
			for (i = 0; i < num_pops; i++)
			{
				last_p = pops[i].p;
				//drift
				binomial_distribution<int> distribution((2 * pops[i].pop_size), pops[i].p);
				pops[i].sampled_p = distribution(generator);
				pops[i].p = (double)pops[i].sampled_p / (2 * (double)pops[i].pop_size);
				//Island model determines p
				//pt = pt-1(1-m)+pbarm
				pops[i].p = pops[i].p*(1 - migration_rate) + (pbar*migration_rate);
				// don't allow wild changes in allele frequencies of 0.8 or higher
				if (pops[i].p - last_p > 0.8 || pops[i].p - last_p < -0.8)
					pops[i].p = last_p;
				if (pops[i].p < 0)
					pops[i].p = 0;
				if (pops[i].p > 1)
					pops[i].p = 1;
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
			ht = 1 - ((total_pop.p*total_pop.p) + (total_pop.q*total_pop.q));
			if (ht > 0)
				fst = (ht - total_pop.het) / ht;
			else
				fst = 0;
			fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
			//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
			//Fst=a/(a+b+c)
			//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
			//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
			//c=0.25hbar
			//nc=(rnbar-sum(n2/rnbar))/(r-1)
			nc = r*nbar;
			s2 = varp = 0;
			for (i = 0; i < num_pops; i++)
			{
				s2 = s2 + pops[i].pop_size*(pops[i].p - total_pop.p)*(pops[i].p - total_pop.p) / ((r - 1)*nbar);
				nc = nc - ((pops[i].pop_size*pops[i].pop_size) / (r*nbar));
				varp = varp + (total_pop.p - pops[i].p)*(total_pop.p - pops[i].p);
			}
			varp = varp / num_pops;
			wc_fst8 = varp / (pbar*(1 - pbar));
			nc = nc / (r - 1);
			pq = (total_pop.p*(1 - total_pop.p));
			rs = ((r - 1) / r)*s2;
			a = (nbar / nc)*(s2 - (1 / (nbar - 1))*(pq - rs - (0.25*total_pop.het)));
			b = (nbar / (nbar - 1))*(pq - rs - ((2 * nbar - 1) / (4 * nbar))*total_pop.het);
			c = 0.5*total_pop.het;
			if ((a + b + c)>0)
				wc_fst = a / (a + b + c);
			else
				wc_fst = 0;

			if (t == 0)
			{
				allele_freqs << "Time\tNumPops\tPbar\tAvgP\tAvgH\tWrightsFst\tWCFst\tWCFstSimple\tNbar\tnc\ts2\tr";
				for (i = 0; i < num_pops; i++)
					allele_freqs << "\tPop" << i << "P\tPop" << i << "Het";
				allele_freqs << '\n' << tt << '\t' << num_pops << '\t' << total_pop.p << '\t' << total_pop.het << '\t'
					<< fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << nbar << '\t' << nc << '\t' << s2 << '\t' << r;
				for (i = 0; i < num_pops; i++)
					allele_freqs << '\t' << pops[i].p << '\t' << pops[i].het;
			}
			
		}//gens
		output << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;
		
		//sample populations
		total_pop.het = total_pop.p = total_pop.q = 0;
		for (i = 0; i < num_pops; i++)
		{
			for (int ii = 0; ii < n_samp_per_pop; ii++)
			{
				sampled_pops.push_back(population()); 
				sampled_pops[samp_index].p = 0;
				for (tt = 0; tt < num_sampled_inds; tt++)
				{
					uniform_real_distribution <double> unidist(0, 1);
					allele1 = unidist(generator);
					allele2 = unidist(generator);
					if (allele1 <= pops[i].p)
						sampled_pops[samp_index].p++;
					if (allele2 <= pops[i].p)
						sampled_pops[samp_index].p++;
				}
				sampled_pops[samp_index].p = sampled_pops[samp_index].p / (2 * num_sampled_inds);
				sampled_pops[samp_index].q = 1 - sampled_pops[samp_index].p;
				sampled_pops[samp_index].het = 2 * sampled_pops[samp_index].p * sampled_pops[samp_index].q;
				total_pop.p = total_pop.p + sampled_pops[samp_index].p;
				total_pop.q = total_pop.q + sampled_pops[samp_index].q;
				total_pop.het = total_pop.het + sampled_pops[samp_index].het;
				samp_index++;
			}
		}
		total_pop.p = total_pop.p / num_sampled;
		total_pop.q = total_pop.q / num_sampled;
		total_pop.het = total_pop.het / num_sampled;
		ht = 1 - ((total_pop.p*total_pop.p) + (total_pop.q*total_pop.q));
		if (ht > 0)
			fst = (ht - total_pop.het) / ht;
		else
			fst = 0;
		fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
		//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
		//Fst=a/(a+b+c)
		//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
		//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
		//c=0.25hbar
		//nc=(rnbar-sum(n2/rnbar))/(r-1)
		nc = r*nbar;
		s2 = varp = 0;
		for (i = 0; i < num_pops; i++)
		{
			s2 = s2 + pops[i].pop_size*(pops[i].p - total_pop.p)*(pops[i].p - total_pop.p) / ((r - 1)*nbar);
			nc = nc - ((pops[i].pop_size*pops[i].pop_size) / (r*nbar));
			varp = varp + (total_pop.p - pops[i].p)*(total_pop.p - pops[i].p);
		}
		varp = varp / num_pops;
		wc_fst8 = varp / (pbar*(1 - pbar));
		nc = nc / (r - 1);
		pq = (total_pop.p*(1 - total_pop.p));
		rs = ((r - 1) / r)*s2;
		a = (nbar / nc)*(s2 - (1 / (nbar - 1))*(pq - rs - (0.25*total_pop.het)));
		b = (nbar / (nbar - 1))*(pq - rs - ((2 * nbar - 1) / (4 * nbar))*total_pop.het);
		c = 0.5*total_pop.het;
		if ((a + b + c)>0)
			wc_fst = a / (a + b + c);
		else
			wc_fst = 0;
		sample_out << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;
	}//reps
	output.close();
	allele_freqs.close();
	sample_out.close();
	
	
	

	if (interactivemode)
	{
		cout << "\nDone! Input integer to quit.\n";
		cin >> end;
		return 0;
	}
	else
	{
		cout << "\n";
		return 0;
	}
}