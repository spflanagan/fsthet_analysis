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
#include "random_numbers.h"

using namespace std;

class population
{
public:
	int pop_size, sampled_p;
	double p, q, het;

	population()
	{
		pop_size = int();
		sampled_p = int();
		p = double();
		q = double();
		het = double();
	}
};

class sampled_inds
{
public:
	vector<int> allele1;
	vector<int> allele2;

	sampled_inds()
	{
		allele1 = vector<int>();
		allele2 = vector<int>();
	}
};

class sampled_pop
{
public:
	vector<sampled_inds> inds;

	sampled_pop()
	{
		inds = vector<sampled_inds>();
	}
};

int main(int argc, char* argv[])
{
	int end, time, i;
	int t, tt, pop_size, sig_count, num_sig;
	int num_pops, num_reps, num_per_p, num_gens, num_sampled, num_sampled_inds;
	double migration_rate, pbar, qbar, fst, ht, last_p, Nm;
	double nc, s2, a, b, c, wc_fst, pq, rs, wc_fst8, varp;
	double w11, w12, w22, s, h;
	int r, nbar, pbar_index;
	bool interactivemode = false;
	bool random_sample = true;
	bool overdominance = false;
	bool directional_selection = false;
	bool sig = false;
	random_device rd;
	default_random_engine generator(rd());
	vector<population> pops;
	population total_pop;
	ofstream allele_freqs, output, sample_out, genepop_out,sig_out, deltaq;
	string base_name, query, tempstring1, tempstring2,sig_name,deltaq_name;
	stringstream output_name, sample_out_name, genepop_out_name;
	vector <bool> pos;

	base_name = "default";	
	s = 0.005;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nnumerical_analysis:\n";
			cout << "Runs numerical analysis to generate Fst-Het patterns\n";
			cout << "-o:\tbase file name (include path). Example: N1000_s10_\n";
			cout << "-n:\tNm\n";
			cout << "-d:\tnumber of demes\n";
			cout << "-r:\trandom sample? In interactive mode, use Y to turn on and N to turn off. In command-line mode, use 1 to turn on.\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-v:\tOverdominance? Follow -v with the selection coefficient (s)\n";
			cout << "-ds:\tDirectional Selection?  Follow -ds with the selection coefficient (s)\n";
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
			cout << "-n:\tNm\n";
			cout << "-d:\tnumber of subpopulations/demes\n";
			cout << "-r:\trandom sample? In interactive mode, use Y to turn on and N to turn off. In command-line mode, use 1 to turn on.\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-v:\tOverdominance? Follow -v with the selection coefficient (s)\n";
			cout << "-ds:\tDirectional Selection?  Follow -ds with the selection coefficient (s)\n";
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
			Nm = atof(tempstring2.c_str());
		if (tempstring1 == "-d")
			num_pops = atoi(tempstring2.c_str());
		if (tempstring1 == "-r")
		{
			if (tempstring2 == "1")
				random_sample = true;
			else
				random_sample = false;
		}
		if (tempstring1 == "-s")
			num_sampled = atoi(tempstring2.c_str());
		if (tempstring1 == "-v")
		{
			overdominance = true;
			s = atof(tempstring2.c_str());
		}
		else
		{
			overdominance = false;
			s = 0;
		}
		if (tempstring1 == "-ds")
		{
			directional_selection = true;
			s = atof(tempstring2.c_str());
		}
		else
		{
			directional_selection = false;
			s = 0;
		}
	}

	if (interactivemode)
	{
		cout << "\nProvide base file name (include path). Example: N1000_s10_\n";
		cin >> base_name;
		cout << "\nProvide Nm\n";
		cin >> Nm;
		cout << "\nProvide the number of subpopulations/demes\n";
		cin >> num_pops;
		cout << "\nProvide the number of populations to sample\n";
		cin >> num_sampled;
		cout << "\nTurn on random sampling? Y or N\n";
		cin >> tempstring2;
		if (tempstring2 == "Y" || tempstring2 == "y")
			random_sample = true;
		else
			random_sample = false;
		cout << "\nInclude Overdominance? Y or N\n";
		cin >> tempstring2;
		if (tempstring2 == "Y" || tempstring2 == "y")
		{
			overdominance = true;
			cout << "\nWhat is the selection coefficient (s)? Suggested is 0.005.\n";
			cin >> s;
		}
		else
			overdominance = false;
		cout << "\nDirectional Selection? Y or N\n";
		cin >> tempstring2;
		if (tempstring2 == "Y" || tempstring2 == "y")
		{
			directional_selection = true;
			cout << "\nWhat is the selection coefficient (s)? Suggested is 0.005.\n";
			cin >> s;
		}
		else
			directional_selection = false;
	}
	
	pop_size = 1000;
	total_pop.pop_size = pop_size * num_pops;

	num_reps = 2000;
	num_sig = 10;
	sig_count = 0;
	num_gens = 5000;
	num_sampled_inds = 50;
	migration_rate = Nm / pop_size;
	cout << "\nRunning with " << num_pops << " pops, " << " sampling " << num_sampled << " of them, with migration rate "
		<< migration_rate << " and a total population size of " << total_pop.pop_size << ".\nRunning " << num_reps 
		<< " number of reps with " << num_gens << " generations.\n";
	if (directional_selection)
		cout << "\nDirectional Selection is imposed on " << num_sig << " loci.\n";
	if (overdominance)
		cout << "\nOverdominance is imposed on " << num_sig << " loci.\n";
	h = 0.5;
	
	for (i = 0; i < num_pops; i++)
	{
		pops.push_back(population());
		if (directional_selection)
			pos.push_back(false);
	}
	if (overdominance || directional_selection)
	{
		sig_name = base_name + "sigloci.txt";
		sig_out.open(sig_name);
		/*deltaq_name = base_name + "delatq.txt";
		deltaq.open(deltaq_name);
		deltaq << "Pop\tGen\tLocus\tObsDeltaq\tDeltaq\tq\tQbar";*/
	}
	if (overdominance)
	{
		w12 = 1;
		w11 = w22 = 1 - s;
	}
	sample_out_name << base_name << "sampledpops.txt";
	sample_out.open(sample_out_name.str());
	sample_out << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	output_name << base_name << "output.txt";
	output.open(output_name.str());
	output << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	stringstream allele_freqs_name;
	/*allele_freqs_name << base_name << "freqs.txt";
	allele_freqs.open(allele_freqs_name.str());
	*/
	vector<population> sampled_pops;
	vector<sampled_pop> samp_pops_genepop;
	double allele1, allele2, total_p, total_het;
	int n_samp_per_pop = (double)num_sampled / (double)num_pops;
	int samp_index = 0;
	for (i = 0; i < num_sampled; i++)
	{
		samp_pops_genepop.push_back(sampled_pop());
		sampled_pops.push_back(population());
		for (int ii = 0; ii < num_sampled_inds; ii++)
		{
			samp_pops_genepop[i].inds.push_back(sampled_inds());
			for (t = 0; t < num_reps; t++)
			{
				samp_pops_genepop[i].inds[ii].allele1.push_back(int());
				samp_pops_genepop[i].inds[ii].allele2.push_back(int());
			}
		}

	}
	for (t = 0; t < num_reps; t++)//each one is a locus
	{
		if (overdominance || directional_selection)
		{
			if (genrand() < 0.5 && sig_count < num_sig)
			{
				if (sig_count == 0)
					sig_out << "loc" << t;
				else
					sig_out << "\nloc" << t;
				//initialize 
				if (directional_selection)
				{
					for (i = 0; i < num_pops; i++)
					{
						if (genrand() < 0.5)
							pos[i] = true;
						else
							pos[i] = false;
						sig_out << '\t' << pos[i];
					}
				}
				sig = true;
				sig_count++;
			}
			else
				sig = false;
		}
		else
			sig = false;
		//set migtrant p-value
		uniform_real_distribution <double> unidist(0.05, 0.95);
		pbar = unidist(generator);
		qbar = 1 - pbar;
		
		if (sig)
		{
			double a,b,c,fourac,qhat1,qhat2, qhat;
			a = (-1*s) / 2;
			b = (s / 2) - migration_rate;
			c = migration_rate*qbar;
			fourac = 4 * a*c;
			qhat1 = ((-1*b) + sqrt(b*b - fourac) )/ (2*a);
			qhat2 = ((-1 * b) -sqrt(b*b - fourac)) / (2 * a);;
			if (qhat1 <= 1 && qhat1 >= 0)
				qhat = qhat1;
			if (qhat2 <= 1 && qhat2 >= 0)
				qhat = qhat2;
			sig_out << '\t' << qbar << '\t' << qhat;
		}

		vector<double>lastps;
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
			double mean_sig_q, mean_null_q, null_count, sig_count;
			mean_sig_q = mean_null_q = null_count = sig_count = 0;
			if (tt == (num_gens - 1) && sig)
				sig_out << "\nloc" << t;
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
				// don't allow wild changes in allele frequencies of 0.8 or higher due todrift
				if (pops[i].p - last_p > 0.8 || pops[i].p - last_p < -0.8)
					pops[i].p = last_p;
				pops[i].q = 1 - pops[i].p;
				//add selection
				if (sig)
				{
					//selection.
					//cout << "\nSelection...p: " << pops[i].p << ", q: " << pops[i].q;
					if (pos[i] == true)
					{
						w11 = 1;
						w12 = 1 - (h*s);
						w22 = 1 - s;

						double qstart = 1 - last_p;
						double wbar = (pops[i].p*pops[i].p*w11) + (2 * pops[i].p*pops[i].q*w12) + (pops[i].q*pops[i].q*w22);
						pops[i].p = ((pops[i].p*pops[i].p*w11) + (pops[i].p*pops[i].q*w12)) / wbar;
						pops[i].q = ((pops[i].q*pops[i].q*w22) + (pops[i].p*pops[i].q*w12)) / wbar;
						double expdq;
						expdq = (s*pops[i].q*pops[i].q)*(1 - pops[i].q)*(1 - pops[i].q + (h*(2 * pops[i].q - 1))) - (migration_rate*(pops[i].q - qbar));
						//deltaq << '\n' << i << '\t' << tt << '\t' << "loc" << t << '\t' << qstart - pops[i].q << '\t' << expdq << '\t' << pops[i].q << '\t' << qbar;							
					}
					if (tt == (num_gens - 1))
					{
						lastps.push_back(last_p);
						sig_out << '\t' << pops[i].q;
						if (pos[i])
						{
							mean_sig_q = mean_sig_q + pops[i].q;
							sig_count++;
						}
						else
						{
							mean_null_q = mean_null_q + pops[i].q;
							null_count++;
						}
					}
				}//end selection
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
			}//end pops
			if (tt == (num_gens - 1) && sig)
				sig_out << '\t' << mean_sig_q / sig_count << '\t' << mean_null_q / null_count;
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
			//fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
			//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
			//Fst=a/(a+b+c)
			//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
			//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
			//c=0.25hbar
			//nc=(rnbar-sum(n2/rnbar))/(r-1)
			r = num_pops;
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

		/*	if (t == 0)
			{
				allele_freqs << "Time\tNumPops\tPbar\tAvgP\tAvgH\tWrightsFst\tWCFst\tWCFstSimple\tNbar\tnc\ts2\tr";
				for (i = 0; i < num_pops; i++)
					allele_freqs << "\tPop" << i << "P\tPop" << i << "Het";
			}
			allele_freqs << '\n' << tt << '\t' << num_pops << '\t' << total_pop.p << '\t' << total_pop.het << '\t'
				<< fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << nbar << '\t' << nc << '\t' << s2 << '\t' << r;
			for (i = 0; i < num_pops; i++)
				allele_freqs << '\t' << pops[i].p << '\t' << pops[i].het;*/
		}//gens
		output << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;
		/*if (sig)
		{
			double qhat;
			for (i = 0; i < num_pops; i++)
			{
				qhat = (s*pops[i].q*pops[i].q)-((migration_rate+s)*pops[i].q)+(migration_rate*qbar);
				deltaq << "\nloc" << t << "\tPop" << i << '\t' << pops[i].q << '\t' << qhat << '\t' << pops[i].p << '\t' << qbar << '\t' << migration_rate << '\t' << s;
			}
		}*/
		
		//sample populations

		
		nbar = 0;
		if (random_sample == false)
		{
			total_pop.het = total_pop.p = total_pop.q = 0;
			samp_index = 0;
			for (i = 0; i < num_pops; i++)
			{
				for (int ii = 0; ii < n_samp_per_pop; ii++)
				{
					sampled_pops[samp_index].pop_size = pops[i].pop_size;
					nbar = nbar + sampled_pops[samp_index].pop_size;
					sampled_pops[samp_index].p = 0;
					for (tt = 0; tt < num_sampled_inds; tt++)
					{
						uniform_real_distribution <double> unidist(0, 1);
						allele1 = unidist(generator);
						allele2 = unidist(generator);
						if (allele1 <= pops[i].p)
						{
							sampled_pops[samp_index].p++;
							samp_pops_genepop[samp_index].inds[tt].allele1[t] = 1;
						}
						else
							samp_pops_genepop[samp_index].inds[tt].allele1[t] = 2;
						if (allele2 <= pops[i].p)
						{
							sampled_pops[samp_index].p++;
							samp_pops_genepop[samp_index].inds[tt].allele2[t] = 1;
						}
						else
							samp_pops_genepop[samp_index].inds[tt].allele2[t] = 2;
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
		}
		if (random_sample == true)
		{
			//establish which pops to sample
			vector<int> samp_indices;
			for (i = 0; i < num_sampled; i++)
			{
				uniform_real_distribution <double> unidist(0, num_pops);
				samp_indices.push_back(int(unidist(generator)));
			}
			total_pop.p = total_pop.q = total_pop.het = 0;
			samp_index = 0;
			for (i = 0; i < num_sampled; i++)
			{
				sampled_pops[i].pop_size = pops[samp_indices[i]].pop_size;
				nbar = nbar + sampled_pops[i].pop_size;
				sampled_pops[i].p = sampled_pops[i].q = sampled_pops[i].het = 0;
				for (tt = 0; tt < num_sampled_inds; tt++)
				{
					uniform_real_distribution <double> unidist(0, 1);
					allele1 = unidist(generator);
					allele2 = unidist(generator);
					if (allele1 <= pops[samp_indices[i]].p)
					{
						sampled_pops[i].p++;
						samp_pops_genepop[i].inds[tt].allele1[t] = 1;
					}
					else
						samp_pops_genepop[i].inds[tt].allele1[t] = 2;
					if (allele2 <= pops[samp_indices[i]].p)
					{
						sampled_pops[i].p++;
						samp_pops_genepop[i].inds[tt].allele2[t] = 1;
					}
					else
						samp_pops_genepop[i].inds[tt].allele2[t] = 2;
				}
				sampled_pops[i].p = sampled_pops[i].p / (2 * num_sampled_inds);
				sampled_pops[i].q = 1 - sampled_pops[i].p;
				sampled_pops[i].het = 2 * sampled_pops[i].p * sampled_pops[i].q;
				total_pop.p = total_pop.p + sampled_pops[i].p;
				total_pop.q = total_pop.q + sampled_pops[i].q;
				total_pop.het = total_pop.het + sampled_pops[i].het;
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
		//fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
		//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
		//Fst=a/(a+b+c)
		//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
		//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
		//c=0.25hbar
		//nc=(rnbar-sum(n2/rnbar))/(r-1)
		r = num_sampled;
		nbar = nbar / num_sampled;
		nc = r*nbar;
		s2 = varp = 0;
		for (i = 0; i < num_sampled; i++)
		{
			s2 = s2 + sampled_pops[i].pop_size*(sampled_pops[i].p - total_pop.p)*(sampled_pops[i].p - total_pop.p) / ((r - 1)*nbar);
			nc = nc - ((sampled_pops[i].pop_size*sampled_pops[i].pop_size) / (r*nbar));
			varp = varp + (total_pop.p - sampled_pops[i].p)*(total_pop.p - sampled_pops[i].p);
		}
		varp = varp / r;
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
//	allele_freqs.close();
	sample_out.close();
	
	if (overdominance || directional_selection)
	{
		sig_out.close();
	//	deltaq.close();
	}

	genepop_out_name << base_name << "genepop";
	genepop_out.open(genepop_out_name.str());
	genepop_out << "Numerical Analysis with Nm=" << Nm << ", N=" << pop_size << ", " << num_pops << " Demes, sampling " << num_sampled << " populations.\nloc0";
	for (t = 1; t < num_reps; t++)
		genepop_out << "\nloc" << t;
	for (t = 0; t < num_sampled; t++)
	{
		genepop_out << "\nPOP"; 
		for (tt = 0; tt < num_sampled_inds; tt++)
		{
			genepop_out << "\nInd" << tt << ",";
			for (i = 0; i < num_reps; i++)
			{
				genepop_out << "\t";
				if (samp_pops_genepop[t].inds[tt].allele1[i] == 1)
					genepop_out << "01";
				else
					genepop_out << "02";
				if (samp_pops_genepop[t].inds[tt].allele2[i] == 1)
					genepop_out << "01";
				else
					genepop_out << "02";
			}
		}
	}
	genepop_out.close();

	
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