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
	int t, tt, pop_size;
	int num_pops, num_reps, num_per_p, num_gens, num_sampled, num_sampled_inds;
	double migration_rate, pbar, qbar, fst, ht, last_p, Nm;
	double nc, s2, a, b, c, wc_fst, pq, rs, wc_fst8, varp;
	double w11, w12, w22;
	int r, nbar, pbar_index;
	bool interactivemode = false;
	bool random_sample = true;
	bool overdominance = false;
	bool directional_selection = false;
	random_device rd;
	default_random_engine generator(rd());
	vector<population> pops;
	population total_pop;
	ofstream allele_freqs, output, sample_out, genepop_out;
	string base_name, query, tempstring1, tempstring2;
	stringstream output_name, sample_out_name, genepop_out_name;

	base_name = "default";	

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
			cout << "-r:\trandom sample? y to turn random sampling on or n to turn random sampling off\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-v:\tOverdominance? Y to turn on or N to turn off\n";
			cout << "-ds:\tDirectional Selection? Y to turn on or N to turn off\n";
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
			cout << "-r:\trandom sample? y to turn random sampling on or n to turn random sampling off\n";
			cout << "-s:\tnumber of populations to sample\n";
			cout << "-v:\tOverdominance? Y to turn on or N to turn off\n";
			cout << "-ds:\tDirectional Selection? Y to turn on or N to turn off\n";
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
			if (tempstring2 == "Y" || tempstring2 == "y")
				random_sample = true;
			else
				random_sample = false;
		}
		if (tempstring1 == "-s")
			num_sampled = atoi(tempstring2.c_str());
		if (tempstring1 == "-v")
		{
			if (tempstring2 == "Y" || tempstring2 == "y")
				overdominance = true;
			else
				overdominance = false;
		}
		if (tempstring1 == "-ds")
		{
			if (tempstring2 == "Y" || tempstring2 == "y")
				directional_selection = true;
			else
				directional_selection = false;
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
			overdominance = true;
		else
			overdominance = false;
		cout << "\nDirectional Selection? Y or N\n";
		cin >> tempstring2;
		if (tempstring2 == "Y" || tempstring2 == "y")
			directional_selection = true;
		else
			directional_selection = false;
	}
	
	pop_size = 1000;
	total_pop.pop_size = pop_size * num_pops;

	num_reps = 2000;
	num_gens = 5000;
	num_sampled_inds = 50;
	migration_rate = Nm / pop_size;
	cout << "\nRunning with " << num_pops << " pops, " << " sampling " << num_sampled << " of them, with migration rate "
		<< migration_rate << " and a total population size of " << total_pop.pop_size << ".\nRunning " << num_reps 
		<< " number of reps with " << num_gens << " generations.\n";
	//initialize variables
	for (i = 0; i < num_pops; i++)
	{
		pops.push_back(population());
	}
	sample_out_name << base_name << "sampledpops.txt";
	sample_out.open(sample_out_name.str());
	sample_out << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	output_name << base_name << "output.txt";
	output.open(output_name.str());
	output << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	stringstream allele_freqs_name;
	allele_freqs_name << base_name << "freqs.txt";
	allele_freqs.open(allele_freqs_name.str());
	
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

			if (t == 0)
			{
				allele_freqs << "Time\tNumPops\tPbar\tAvgP\tAvgH\tWrightsFst\tWCFst\tWCFstSimple\tNbar\tnc\ts2\tr";
				for (i = 0; i < num_pops; i++)
					allele_freqs << "\tPop" << i << "P\tPop" << i << "Het";
			}
			allele_freqs << '\n' << tt << '\t' << num_pops << '\t' << total_pop.p << '\t' << total_pop.het << '\t'
				<< fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << nbar << '\t' << nc << '\t' << s2 << '\t' << r;
			for (i = 0; i < num_pops; i++)
				allele_freqs << '\t' << pops[i].p << '\t' << pops[i].het;
		}//gens
		output << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;
		
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
	allele_freqs.close();
	sample_out.close();

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