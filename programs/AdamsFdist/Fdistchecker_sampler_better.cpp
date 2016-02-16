#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>

using namespace std;

class population
{
public:
	double* freq_p;
	double* h_s; 
};

class sampled_alleles
{
public:
	vector< vector<char> > allele_list;
};

int main()
{
	unsigned seed = static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count());
	default_random_engine generator(seed);

	int numberloci;
	double m;
	int numberpops;
	double* migrant_pool_p;
	int pop_size;
	int loci_loop;
	int pop_loop;
	int gens_loop;
	int numbergenerations;
	double mean_p;
	vector<double> pop_p;
	vector<double> pop_h;
	
	int i, j, k;
	int number_samples_per_pop, sample_size, samplecounter;
	vector<double> sample_p;
	vector<double> sample_hs;

	bool randomly_sample;
	int numberpops_randomsample;

	bool output_genepop;

	// set parameters

	m = 0.01;
	numberpops = 50;
	numberpops_randomsample = 20;// only used if populations are sampled at random
	
	randomly_sample = true;
	output_genepop = true;
	pop_size = 1000;
	numbergenerations = 5000;
	numberloci = 2000;
	number_samples_per_pop = 2; // only used if every pop is sampled the same number of times
	sample_size = 50;

	population *pop;
	pop = new population[numberpops];
	for (i = 0; i < numberpops; i++)
	{
		pop[i].freq_p = new double[numberloci];
		pop[i].h_s = new double[numberloci];
	}
	migrant_pool_p = new double[numberloci];

	for (i = 0; i < numberpops; i++)
	{
		for (j = 0; j < numberloci; j++)
		{
			pop[i].freq_p[j] = 0.5;
			pop[i].h_s[j] = 0.5;
		}
	}

	// set the migrant pool p for all locus here, which effectively sets the heterozygosity
	uniform_real_distribution<double> unidist(0.05,0.95);

	for (i = 0; i < numberloci; i++)
		migrant_pool_p[i] = unidist(generator);

	for (gens_loop = 0; gens_loop < numbergenerations; gens_loop++)
	{
		double dHs = 0;
		mean_p = 0;
		for (pop_loop = 0; pop_loop < numberpops; pop_loop++)
		{
			for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
			{
				double p_before = pop[pop_loop].freq_p[loci_loop];
				binomial_distribution<int> distribution(2*pop_size,pop[pop_loop].freq_p[loci_loop]);
			// drift
				int iAlleleCount = distribution(generator);
				double dAlleleCount = static_cast<double>(iAlleleCount);
				double dTotAlleles = static_cast<double>(2*pop_size);
				pop[pop_loop].freq_p[loci_loop] = dAlleleCount/dTotAlleles;

				// don't allow wild changes in allele frequencies of 0.8 or higher
				if (pop[pop_loop].freq_p[loci_loop]-p_before > 0.8 || pop[pop_loop].freq_p[loci_loop]-p_before < -0.8)
					pop[pop_loop].freq_p[loci_loop] = p_before;
				// make sure p is between 0 and 1
				if (pop[pop_loop].freq_p[loci_loop] > 1)
					pop[pop_loop].freq_p[loci_loop] = 1;
				if (pop[pop_loop].freq_p[loci_loop] < 0)
					pop[pop_loop].freq_p[loci_loop] = 0;

				// migration
				pop[pop_loop].freq_p[loci_loop] = (pop[pop_loop].freq_p[loci_loop] * (1.0-m)) + (migrant_pool_p[loci_loop] * m);

			} // loci loop
		} // pop_loop		
	} // gens_loop

	// Sample from each population
	
	double dsamplesize = static_cast<int>(sample_size);
	uniform_real_distribution<double> unidist1(0.0,1.0);

	sampled_alleles* allele_sample;
	population* sample;

	int total_number_of_samples;
	if (!randomly_sample)
		total_number_of_samples = number_samples_per_pop*numberpops;
	else
		total_number_of_samples = numberpops_randomsample;

	sample = new population[total_number_of_samples];
	for (i = 0; i < total_number_of_samples; i++)
	{
		sample[i].freq_p = new double[numberloci];
		sample[i].h_s = new double[numberloci];
	}

	allele_sample = new sampled_alleles[total_number_of_samples];
	vector<char> temp_allele_list;

	if (!randomly_sample)  // Here at least one sample is taken from every deme
	{
		samplecounter = 0;
		for (i = 0; i < numberpops; i++)
		{
			for (j = 0; j < number_samples_per_pop; j++)
			{
				for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
				{
					temp_allele_list.clear();
					sample[samplecounter].freq_p[loci_loop] = 0;
					for (k = 0; k < sample_size*2; k++)
					{
						if (unidist1(generator) < pop[i].freq_p[loci_loop])
						{
							sample[samplecounter].freq_p[loci_loop] = sample[samplecounter].freq_p[loci_loop] + 1;
							temp_allele_list.push_back('A');
						}
						else
						{
							temp_allele_list.push_back('a');
						}
					} //k
					sample[samplecounter].freq_p[loci_loop] = sample[samplecounter].freq_p[loci_loop]/(dsamplesize*2);
					allele_sample[samplecounter].allele_list.push_back(temp_allele_list);
				} // loci_loop
				samplecounter++;
			} // j
		} // i
	}
	else  // Here a number of populations are sampled at random from the demes
	{
		int currentpop;
		double dcurrentpop;
		samplecounter = 0;
		for (i = 0; i < numberpops_randomsample; i++)
		{
			dcurrentpop = unidist1(generator) * numberpops;
			currentpop = static_cast<int>(floor(dcurrentpop));

			for(loci_loop = 0; loci_loop < numberloci; loci_loop++)
			{
				temp_allele_list.clear();
				sample[samplecounter].freq_p[loci_loop] = 0;
				for (k = 0; k < sample_size*2; k++)
				{
					if (unidist1(generator) < pop[currentpop].freq_p[loci_loop])
					{
						sample[samplecounter].freq_p[loci_loop] = sample[samplecounter].freq_p[loci_loop] + 1;
						temp_allele_list.push_back('A');
					}
					else
					{
						temp_allele_list.push_back('a');
					}
				} //k
				sample[samplecounter].freq_p[loci_loop] = sample[samplecounter].freq_p[loci_loop]/(dsamplesize*2);	
				allele_sample[samplecounter].allele_list.push_back(temp_allele_list);
			}
			samplecounter++;
		}
	}

	// Now the population loci all have allele frequencies
	// The sampled loci also have allele frequencies
	// Calculate Fst, Hs, Ht

	double dnumsamples = static_cast<double>(samplecounter);
	double dnumpops = static_cast<double>(numberpops);

	// for every locus, we will need a mean Hs and a mean p in the sample and the population

	double *sample_mean_hs = new double[numberloci];
	double *sample_mean_p = new double[numberloci];
	double *sample_overall_ht = new double[numberloci];
	double *sample_overall_fst = new double[numberloci];
	double *pop_mean_hs = new double[numberloci];
	double *pop_mean_p = new double[numberloci];
	double *pop_overall_ht = new double[numberloci];
	double *pop_overall_fst = new double[numberloci];

	for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
	{
		sample_mean_p[loci_loop] = 0;
		sample_mean_hs[loci_loop] = 0;
		pop_mean_p[loci_loop] = 0;
		pop_mean_hs[loci_loop] = 0;

		for (i = 0; i < numberpops; i++)
		{
			pop[i].h_s[loci_loop] = 2*pop[i].freq_p[loci_loop]*(1-pop[i].freq_p[loci_loop]);
			pop_mean_p[loci_loop] = pop_mean_p[loci_loop] + pop[i].freq_p[loci_loop];
			pop_mean_hs[loci_loop] = pop_mean_hs[loci_loop] + pop[i].h_s[loci_loop];
		}
		pop_mean_p[loci_loop] = pop_mean_p[loci_loop]/dnumpops;
		pop_mean_hs[loci_loop] = pop_mean_hs[loci_loop]/dnumpops;
		pop_overall_ht[loci_loop] = 2*pop_mean_p[loci_loop]*(1-pop_mean_p[loci_loop]);
		if (pop_overall_ht[loci_loop] > 0)
			pop_overall_fst[loci_loop] = 1-(pop_mean_hs[loci_loop]/pop_overall_ht[loci_loop]);
		else
			pop_overall_fst[loci_loop] = 0;

		for (i = 0; i < samplecounter; i++)
		{
			sample[i].h_s[loci_loop] = 2*sample[i].freq_p[loci_loop]*(1-sample[i].freq_p[loci_loop]);
			sample_mean_p[loci_loop] = sample_mean_p[loci_loop] + sample[i].freq_p[loci_loop];
			sample_mean_hs[loci_loop] = sample_mean_hs[loci_loop] + sample[i].h_s[loci_loop];
		}
		sample_mean_p[loci_loop] = sample_mean_p[loci_loop]/dnumsamples;
		sample_mean_hs[loci_loop] = sample_mean_hs[loci_loop]/dnumsamples;
		sample_overall_ht[loci_loop] = 2*sample_mean_p[loci_loop]*(1-sample_mean_p[loci_loop]);
		if (sample_overall_ht[loci_loop] > 0)
			sample_overall_fst[loci_loop] = 1-(sample_mean_hs[loci_loop]/sample_overall_ht[loci_loop]);
		else
			sample_overall_fst[loci_loop] = 0;
	} // loci_loop


	// Calculate variance-based Fst
	double* sample_var_fst = new double[numberloci];
	double varianceinp;
	for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
	{
		varianceinp = 0;
		for (i = 0; i < samplecounter; i++)
		{
			varianceinp = varianceinp + (sample[i].freq_p[loci_loop]-sample_mean_p[loci_loop])*(sample[i].freq_p[loci_loop]-sample_mean_p[loci_loop]);
		}
		varianceinp = varianceinp/(dnumsamples);
		if (sample_mean_p[loci_loop] > 0 && sample_mean_p[loci_loop] < 1)
			sample_var_fst[loci_loop] = varianceinp/((sample_mean_p[loci_loop])*(1-sample_mean_p[loci_loop]))-1/(2*dsamplesize);
		else
			sample_var_fst[loci_loop] = 0;
	}

	// Calculate the Weir and Cockerham Fst (Genetic Data Analysis, p. 147-148)
	// I'm using the haploid estimator because we assume everything's in HWE within pops.

	double* sample_theta = new double[numberloci];
	double theta;
	double T_1, T_2;
	double dnbar, dn_c, d_r, d_s2a, dn_i;

	for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
	{
		dn_i = static_cast<double>(sample_size);
		d_r = dnumsamples;
		dnbar = dn_i; // the mean n equals the n of population i, because all samples are the same size

		// Calculate the variance in allele frequencies
		d_s2a = 0;
		for (i = 0; i < samplecounter; i++)
		{
			d_s2a = d_s2a + dn_i*(sample[i].freq_p[loci_loop]-sample_mean_p[loci_loop])*(sample[i].freq_p[loci_loop]-sample_mean_p[loci_loop]);
		}
		d_s2a = d_s2a/((d_r-1)*dnbar);

		// Calculate the strange n(c) term
		double sum_of_ni2;
		double sum_of_ni;
		sum_of_ni = 0;
		sum_of_ni2 = 0;
		for (i = 0; i < samplecounter; i++)
		{
			sum_of_ni = sum_of_ni + dn_i;
			sum_of_ni2 = sum_of_ni2 + dn_i*dn_i;
		}
		dn_c = (1/(d_r-1))*(sum_of_ni - sum_of_ni2/sum_of_ni);

		// Calculate theta
		T_1 = d_s2a - (1/(dnbar-1))*(sample_mean_p[loci_loop]*(1-sample_mean_p[loci_loop])-((d_r-1)/d_r)*d_s2a);
		T_2 = ((dn_c-1)/(dnbar-1))*(sample_mean_p[loci_loop]*(1-sample_mean_p[loci_loop]))+(d_s2a/d_r)*(1 + (d_r-1)*(dnbar-dn_c)/(dnbar-1));
		if (T_2 != 0)
			theta = T_1/T_2;
		else
			theta = 0;

		sample_theta[loci_loop] = theta;
	}


	ofstream result_file;
	ofstream genepop_file;

	result_file.open("Fdistcheckerresults.xls");

	result_file << "Loc\tm\tN\tN_pops\tsamples\ts_size\tm_p\tHs\tHt\tFst\ttheta\tvarFst\tpopHt\tpopFst";
	for (loci_loop = 0; loci_loop < numberloci; loci_loop++)
	{
		result_file << "\n" << loci_loop + 1 << "\t" << m << "\t" << pop_size << "\t";
		result_file << numberpops << "\t";
		result_file << samplecounter << "\t";
		result_file << sample_size << "\t" << migrant_pool_p[loci_loop] << "\t" << sample_mean_hs[loci_loop] << "\t";
		result_file << sample_overall_ht[loci_loop] << "\t" << sample_overall_fst[loci_loop] << "\t";
		result_file << sample_theta[loci_loop] << "\t" << sample_var_fst[loci_loop] << "\t";
		result_file << pop_overall_ht[loci_loop] << "\t" << pop_overall_fst[loci_loop];
	}
	result_file.close();
	
	if (output_genepop)
	{
		genepop_file.open("genepopdata.txt");
		genepop_file << "TitleLine";
		for (i = 0; i < numberloci; i++)
		{
			genepop_file << "\nLocus" << i+1;
		}

		for (i = 0; i < samplecounter; i++)
		{
			genepop_file << "\npop";
			for (j = 0; j < sample_size; j++)
			{
				genepop_file << "\nind" << j+1 << ",";
				for (k = 0; k < numberloci; k++)
				{
					genepop_file << " ";
					if (allele_sample[i].allele_list[k][j*2] == 'A')
						genepop_file << "13";
					else
						genepop_file << "15";
					if (allele_sample[i].allele_list[k][j*2+1] == 'A')
						genepop_file << "13";
					else
						genepop_file << "15";
				}
			}
		}
		genepop_file.close();
	}	


	for (i = 0; i < numberpops; i++)
	{
		delete[] pop[i].freq_p;
		delete[] pop[i].h_s;
	}
	delete[] pop;

	delete[] migrant_pool_p;

		
	for (i = 0; i < samplecounter; i++)
	{
		delete[] sample[i].freq_p;
		delete[] sample[i].h_s;
	}
	delete[] sample;

	delete[] sample_mean_hs;
	delete[] sample_mean_p;
	delete[] sample_overall_ht;
	delete[] sample_overall_fst;
	delete[] pop_mean_hs;
	delete[] pop_mean_p;
	delete[] pop_overall_ht;
	delete[] pop_overall_fst;
	delete[] sample_var_fst;
	delete[] sample_theta;
	delete[] allele_sample;

	return 0;

}