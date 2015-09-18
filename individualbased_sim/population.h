#pragma once
#include "random_numbers.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <vector>
#include <array>
#include <sstream>

using namespace std;

class parameters
{
public:
	int num_reps, carrying_capacity;
	int max_fecundity, max_encounters, max_progeny;
	int total_num_loci, num_poly_loci;
	int sample_size, generations, num_pops;
	int num_chrom, num_markers, num_qtls, num_alleles;
	double recomb_rate, mut_rate, mut_var, env_var;
	double GaussianPrefMean, gauss_var, migration_rate, viability_sel;

	void set_parameters()
	{
		//Initialize Parameters
		num_reps = 1;
		carrying_capacity = 1000;
		max_fecundity = 4;
		max_progeny = carrying_capacity*max_fecundity;
		num_chrom = 2;
		num_markers = 100;
		num_qtls = 2;
		GaussianPrefMean = 4.0;
		num_alleles = 4;
		max_encounters = 50;
		recomb_rate = 0.2;
		mut_rate = 0;
		mut_var = 0;
		total_num_loci = num_markers*num_chrom;
		env_var = 0;
		sgenrand(time(0));
		sample_size = 50;//set sample size here
		generations = 2000;
		num_pops = 5;
		gauss_var = 0; //if 0, it's random mating
		viability_sel = 0;
		migration_rate = 0.1;
	}

};

class reference_locus_info
{
public:
	vector <char> snps;
	vector <double> snp_freq;

	reference_locus_info()
	{
		snps = vector <char>();
		snp_freq = vector<double>();
	}
};

class chromosome
{
public:
	vector<int> loci;
	vector<double> allelic_effects;
};

class individual
{
public:
	vector<chromosome> maternal;
	vector<chromosome> paternal;
	double phenotype, genotype;
	bool alive;//true=alive,false=dead
	bool female;//true=female, false=male
	int mate_found;

	individual()//class constructor
	{
		maternal = vector < chromosome >() ;
		paternal = vector<chromosome>();
		phenotype = double();
		genotype = double();
		alive = bool();
		female = bool();
		mate_found = bool();
	}
};//end Individual

class population
{
public:
	vector<individual> adult, progeny;
	vector<vector<int>> qtls_locations;
	int population_size, num_males, num_females, num_progeny, num_not_mated;
	double mean_mal_trait, mean_fem_trait;
	bool extinct;

public:
	void initialize(vector <int> &qtls, parameters &p)
	{
		extinct = false;
		int j,jj,jjj, count,q, qtl_counter;
		//set up qtl tracker
		count = q = qtl_counter = 0;
		for (j = 0; j < p.num_chrom; j++)
		{
			qtl_counter = 0;
			qtls_locations.push_back(vector<int>());
			for (jj = 0; jj < p.num_markers; jj++)
			{
				if (count == qtls[q])
				{
					qtls_locations[j].push_back(qtl_counter);
					qtl_counter++;
					q++;
				}
				else
					qtls_locations[j].push_back(-1);
				count++;
			}
		}
		//set up genotypes
		num_females = num_males = num_progeny = 0;
		for (j = 0; j < p.carrying_capacity; j++)
		{
			adult.push_back(individual());
			if (genrand() < 0.5)
			{
				adult[j].female = true;
				num_females++;
			}
			else
			{
				adult[j].female = false;
				num_males++;
			}
			adult[j].alive = true;
			for (jj = 0; jj < p.num_chrom; jj++)
			{
				adult[j].maternal.push_back(chromosome());
				adult[j].paternal.push_back(chromosome());
				for (jjj = 0; jjj < p.num_markers; jjj++)
				{
					adult[j].maternal[jj].loci.push_back(randnum(p.num_alleles));
					adult[j].paternal[jj].loci.push_back(randnum(p.num_alleles));
				}
			}
		}
		//assign allelic standard deviations
		vector<double> temp_allele;
		double allelic_sd = 0.5;
		for (j = 0; j < p.num_alleles; j++)
			temp_allele.push_back(randnorm(0, allelic_sd));
		for (j = 0; j < p.carrying_capacity; j++)
		{
			adult[j].phenotype = 0;
			adult[j].genotype = 0;
			for (jj = 0; jj < p.num_chrom; jj++)
			{
				//Assign allelic effects
				for (jjj = 0; jjj < qtls_locations[jj].size(); jjj++)
				{
					adult[j].maternal[jj].allelic_effects.push_back(temp_allele[j%p.num_alleles]);
					adult[j].paternal[jj].allelic_effects.push_back(temp_allele[j%p.num_alleles]);
					adult[j].genotype = adult[j].genotype +
						adult[j].maternal[jj].allelic_effects[jjj] +
						adult[j].paternal[jj].allelic_effects[jjj];
				}
				adult[j].phenotype = adult[j].genotype + randnorm(0, sqrt(p.env_var));
			}
			if (adult[j].female)
				mean_fem_trait = mean_fem_trait + adult[j].phenotype;
			else
				mean_mal_trait = mean_mal_trait + adult[j].phenotype;
		}
		population_size = p.carrying_capacity;
		mean_mal_trait = mean_mal_trait / num_males;
		mean_fem_trait = mean_fem_trait / num_females;
		
	}//initialize

	void RecombineChromosome(chromosome &RecombinedChr, individual &Parent, int WhichChromosome, double ExpectedRecombEvents, parameters &p)
	{
		int RCi, RCj;
		int NumberRecombEvents = 0;
		int SegmentStart[22], SegmentEnd[22];
		int BreakPoint[20];

		if (ExpectedRecombEvents < 6)
			NumberRecombEvents = poissonrand(ExpectedRecombEvents);
		if (ExpectedRecombEvents > 6)
			NumberRecombEvents = positiveroundnorm(ExpectedRecombEvents, sqrt(ExpectedRecombEvents));
		for (RCi = 0; RCi < 20; RCi++)
			BreakPoint[RCi] = p.num_markers + 1;
		bool SegmentMaternal[22];
		int NumberSegments;
		bool StartMaternal;
		if (NumberRecombEvents > 20)
			NumberRecombEvents = 20;

		if (NumberRecombEvents > 0)
		{
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
				BreakPoint[RCi] = randnum(p.num_markers);
			//sort breakpoints
			sort(begin(BreakPoint), end(BreakPoint));

			//first segment maternal or paternal?
			if (genrand() < 0.5)
				StartMaternal = true;
			else
				StartMaternal = false;

			NumberSegments = 1;
			SegmentStart[0] = 0;
			SegmentMaternal[0] = StartMaternal;
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
			{
				SegmentEnd[RCi] = BreakPoint[RCi];
				SegmentStart[RCi + 1] = BreakPoint[RCi];
				if (SegmentMaternal[RCi])
					SegmentMaternal[RCi + 1] = false;
				else
					SegmentMaternal[RCi + 1] = true;
				NumberSegments++;
			}//end RCi
			SegmentEnd[RCi] = p.num_markers;

			//now pass allelic info to recombined chromosome
			for (RCi = 0; RCi < NumberSegments; RCi++)
			{
				if (SegmentMaternal[RCi])
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
					{
						RecombinedChr.loci[RCj] = Parent.maternal[WhichChromosome].loci[RCj];
						if (qtls_locations[WhichChromosome][RCj] >= 0)
							RecombinedChr.allelic_effects[qtls_locations[WhichChromosome][RCj]] = Parent.maternal[WhichChromosome].allelic_effects[qtls_locations[WhichChromosome][RCj]];
					}
				}
				else
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
					{
						RecombinedChr.loci[RCj] = Parent.paternal[WhichChromosome].loci[RCj];
						if (qtls_locations[WhichChromosome][RCj] >= 0)
							RecombinedChr.allelic_effects[qtls_locations[WhichChromosome][RCj]] = Parent.paternal[WhichChromosome].allelic_effects[qtls_locations[WhichChromosome][RCj]];
					}
				}
			}//end RCi
		}//end numb recomb events > 0
		else
		{
			//No recombination
			if (genrand() < 0.5)
			{
				for (RCi = 0; RCi < p.num_markers; RCi++)
					RecombinedChr.loci[RCi] = Parent.maternal[WhichChromosome].loci[RCi];
				for (RCi = 0; RCi < qtls_locations[WhichChromosome].size(); RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.maternal[WhichChromosome].allelic_effects[RCi];
			}
			else
			{
				for (RCi = 0; RCi < p.num_markers; RCi++)
					RecombinedChr.loci[RCi] = Parent.paternal[WhichChromosome].loci[RCi];
				for (RCi = 0; RCi < qtls_locations[WhichChromosome].size(); RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.paternal[WhichChromosome].allelic_effects[RCi];
			}
		}//else (no recomb)
	}//end RecombineChromosome

	void AssignProgenyGenotypes(individual &prog, individual &parent, bool maternal, parameters &p)
	{
		int j, jj, jjj;
		if (genrand() > 0.5)//goes to offspring maternal
		{
			for (j = 0; j < p.num_chrom; j++)
			{
				for (jj = 0; jj < p.num_markers; jj++)
				{
					if (maternal)
					{
						prog.maternal[j].loci[jj] = parent.maternal[j].loci[jj];
						for (jjj = 0; jjj < qtls_locations[j].size(); jjj++)
						{
							prog.maternal[j].allelic_effects[jjj] = parent.maternal[j].allelic_effects[jjj];
						}
					}
					else
					{
						prog.maternal[j].loci[jj] = parent.paternal[j].loci[jj];
						for (jjj = 0; jjj < qtls_locations[j].size(); jjj++)
						{
							prog.maternal[j].allelic_effects[jjj] = parent.paternal[j].allelic_effects[jjj];
						}
					}					
				}//end of num_markers
			}
		}
		else
		{
			for (j = 0; j < p.num_chrom; j++)
			{
				for (jj = 0; jj < p.num_markers; jj++)
				{
					if (maternal)
					{
						prog.paternal[j].loci[jj] = parent.maternal[j].loci[jj];
						for (jjj = 0; jjj < qtls_locations[j].size(); jjj++)
						{
							prog.paternal[j].allelic_effects[jjj] = parent.maternal[j].allelic_effects[jjj];
						}
					}
					else
					{
						prog.paternal[j].loci[jj] = parent.paternal[j].loci[jj];
						for (jjj = 0; jjj < qtls_locations[j].size(); jjj++)
						{
							prog.paternal[j].allelic_effects[jjj] = parent.paternal[j].allelic_effects[jjj];
						}
					}
				}//end of num_markers
			}
		}
	}

	void Mating(parameters &p)
	{
		int j, jj, jjj, k, irndnum;
		int num_prog, females, males, fem_prog, mal_prog, num_mated;
		double dubrand;
		double mean_fem_trait, mean_mal_trait, mate_prob;
		bool mate_found, males_present;
		int mal_index, fem_id, mal_id, encounters;
		int counter1, counter2, counter3;
		vector <int> male_list, fathers, mothers, non_mated;
		num_prog = females = males = fem_prog = mal_prog = num_mated = 0;
		males_present = false;
		mean_mal_trait = mean_fem_trait = 0;
		double sd_mal_trait = 0;
		counter1 = counter2 = counter3 = 0;

		UpdateDemographics();
		progeny.resize(0);
		//determine mean male trait
		for (j = 0; j < population_size; j++)
		{
			adult[j].mate_found = 0;
			if (!adult[j].female)
			{
				males_present = true;
				male_list.push_back(j);
				males++;
				mean_mal_trait = mean_mal_trait+ adult[j].phenotype;
			}
		} // end of males
		if (males > 0) {
			mean_mal_trait = mean_mal_trait/ males;
		}
		else {
			mean_mal_trait = 0;
			extinct = true;
		}
		for (j = 0; j < population_size; j++){
			if (!adult[j].female)
				sd_mal_trait = sd_mal_trait + (adult[j].phenotype - mean_mal_trait)*(adult[j].phenotype - mean_mal_trait);
		}
		sd_mal_trait = sqrt(sd_mal_trait/ males);
		for (j = 0; j < population_size; j++)
		{
			mate_found = false;
			if (adult[j].female && males_present)
			{
				fem_id = j;
				mean_fem_trait = mean_fem_trait + adult[j].phenotype;
				females++;
				encounters = 0;
				while (!mate_found && encounters <= p.max_encounters)
				{
					irndnum = randnum(males);
					mal_index = male_list[irndnum];
					if (p.gauss_var > 0)
						mate_prob = exp(-0.5 * (adult[mal_index].phenotype - p.GaussianPrefMean)*
						(adult[mal_index].phenotype - p.GaussianPrefMean) / p.gauss_var);
					else//it's random mating
						mate_prob = 1;
					dubrand = genrand();
					if (dubrand < mate_prob)
					{
						mate_found = true;
						mal_id = mal_index;
						fathers.push_back(mal_id); //counter2
						adult[mal_id].mate_found++;
						counter2++;
						num_mated++;
					}
					encounters++;
				}//while
				if (mate_found)//then they mate!
				{
					mothers.push_back(j); //counter1
					counter1++;
					adult[j].mate_found++;
					if (num_prog >= p.max_progeny)
						num_prog = p.max_progeny - 1;
					//mother is Parent 1, j
					//father is parent 2, mal_id
					for (k = 0; k < p.max_fecundity; k++)
					{
						progeny.push_back(individual());
						for (jj = 0; jj < p.num_chrom; jj++)
						{
							progeny[num_prog].maternal.push_back(chromosome());
							progeny[num_prog].paternal.push_back(chromosome());
							for (jjj = 0; jjj < p.num_markers; jjj++)
							{
								progeny[num_prog].maternal[jj].loci.push_back(int());
								progeny[num_prog].paternal[jj].loci.push_back(int());
							}
							for (jjj = 0; jjj < qtls_locations[jj].size(); jjj++)
							{
								progeny[num_prog].maternal[jj].allelic_effects.push_back(double());
								progeny[num_prog].paternal[jj].allelic_effects.push_back(double());
							}
						}
						AssignProgenyGenotypes(progeny[num_prog], adult[fem_id], true, p);
						AssignProgenyGenotypes(progeny[num_prog], adult[mal_id], false, p);
						progeny[num_prog].alive = true;
						//calculate phenotype in mutation once the genotype is determined
						for (jjj = 0; jjj < p.num_chrom; jjj++)//go through chromosome by chromosome
						{
							RecombineChromosome(progeny[num_prog].maternal[jjj], adult[fem_id], jjj, p.recomb_rate, p);
							RecombineChromosome(progeny[num_prog].paternal[jjj], adult[mal_id], jjj, p.recomb_rate, p);
						}//end of chromosome
						if (genrand() < 0.5)
						{
							progeny[num_prog].female = true;
							fem_prog++;
						}
						else
						{
							progeny[num_prog].female = false;
							mal_prog++;
						}
						num_prog++;
					}//for k
				}
				if (!mate_found)
				{//Keep track of the females that didn't mate
					non_mated.push_back(j);//counter3
					counter3++;
				}
			}//if female
		}//end of mm
		num_progeny = num_prog;
		num_females = females;
		num_males = males;
		num_not_mated = counter3;
		mean_fem_trait = mean_fem_trait / females;
		progeny.resize(num_progeny);
	}//end mating

	void Mutation(parameters &p)
	{
		int j, jj, jjj, irand, irand2, irand3, locus;
		double rnd1, rnd2;
		double ind_mut_rate, mut_sd;
		bool mutated;

		mut_sd = sqrt(p.mut_var);
		ind_mut_rate = p.mut_rate * 2 * p.total_num_loci;

		for (j = 0; j < num_progeny; j++)
		{
			rnd1 = genrand();
			progeny[j].phenotype = 0;
			progeny[j].genotype = 0;
			mutated = false;
			if (rnd1 < ind_mut_rate)
			{
				irand = randnum(p.num_chrom);
				irand2 = randnum(p.num_markers);
				rnd2 = genrand();//to choose maternal or paternal
				if (rnd2 < 0.5)//affects maternal chromosome	
				{
					while (!mutated){
						irand3 = randnum(p.num_alleles);
						if (!progeny[j].maternal[irand].loci[irand2] == irand3)
						{
							progeny[j].maternal[irand].loci[irand2] = irand3;
							if (qtls_locations[irand][irand2] >= 0)
							{
								progeny[j].maternal[irand].allelic_effects[qtls_locations[irand][irand2]] =
									progeny[j].maternal[irand].allelic_effects[qtls_locations[irand][irand2]] + randnorm(0, mut_sd);
							}
							mutated = true;
						}
					}
				}
				else//affects paternal chromosome
				{
					while (!mutated){
						irand3 = randnum(p.num_alleles);
						if (!progeny[j].paternal[irand].loci[irand2] == irand3)
						{
							progeny[j].paternal[irand].loci[irand2] = irand3;
							if (qtls_locations[irand][irand2] >= 0)
							{
								progeny[j].maternal[irand].allelic_effects[qtls_locations[irand][irand2]] =
									progeny[j].maternal[irand].allelic_effects[qtls_locations[irand][irand2]] + randnorm(0, mut_sd);
							}
							mutated = true;
						}
					}
				}
			}//end of if

			for (jj = 0; jj < p.num_chrom; jj++)
			{
				for (jjj = 0; jjj < p.num_markers; jjj++){
					progeny[j].genotype = progeny[j].genotype +
						progeny[j].maternal[jj].allelic_effects[jjj] + progeny[j].paternal[jj].allelic_effects[jjj];
				}
			}//end of jj
			progeny[j].phenotype = progeny[j].genotype + randnorm(0, sqrt(p.env_var));
		}//end of j
	}//mutation

	void SelectionOnPhenotypes(double dSelectionStrength)
	{
		int j, prog_alive, male_count;
		double survival_prob, drnum1, optimum, phen_sd, phen_m, num;
		male_count = prog_alive = 0;
		phen_sd = phen_m = 0;
		//calc mean male phenotype & std dev
		for (j = 0; j < num_progeny; j++){
			if (!progeny[j].female){
				phen_m = phen_m + progeny[j].phenotype;
				male_count++;
			}
		}
		num = male_count;
		phen_m = phen_m / num;
		for (j = 0; j < num_progeny; j++){
			if (!progeny[j].female)
				phen_sd = phen_sd + (progeny[j].phenotype - phen_m)*(progeny[j].phenotype - phen_m);
		}
		phen_sd = sqrt(phen_sd/ num);

		optimum = 0;
		for (j = 0; j < num_progeny; j++) {
			if (progeny[j].female)
			{
				progeny[j].alive = true;
				prog_alive++;
			}
			else//selection only on males
			{
				if (dSelectionStrength > 0)
					survival_prob = exp(-1 * (progeny[j].phenotype - optimum)*(progeny[j].phenotype - optimum)
					/ (2 * dSelectionStrength));
				else
					survival_prob = 1;
				//cout<<dSurvProb<<'\n';
				drnum1 = genrand();
				if (drnum1 < survival_prob)
				{
					progeny[j].alive = true;
					prog_alive++;
				}
				else
					progeny[j].alive = false;
			}
		} // end of i
	}//end selection on phenotypes

	void DensityRegulation(parameters &p)
	{
		int j, jj, jjj, num_adults_chosen;
		double car_cap_rem, prog_left, keep_prob, drnd;
		num_males = num_females = prog_left = 0;
		//count the ones that are still alive
		for (j = 0; j < num_progeny; j++)
		{
			if (progeny[j].alive)
				prog_left++;
		}
		car_cap_rem = p.carrying_capacity;
		num_adults_chosen = 0;
		for (j = 0; j < num_progeny; j++)
		{
			if (progeny[j].alive)
			{
				if (prog_left == 0)
					keep_prob = 0;
				else
					keep_prob = car_cap_rem / prog_left;
				drnd = genrand();
				if (drnd < keep_prob)
				{//then turn it into an adult
					adult[num_adults_chosen].alive = true;
					adult[num_adults_chosen].mate_found = 0;
					for (jj = 0; jj < p.num_chrom; jj++)
					{
						for (jjj = 0; jjj < p.num_markers; jjj++)
						{
							adult[num_adults_chosen].maternal[jj].loci[jjj] = progeny[j].maternal[jj].loci[jjj];
							adult[num_adults_chosen].paternal[jj].loci[jjj] = progeny[j].paternal[jj].loci[jjj];
							if (qtls_locations[jj][jjj] >= 0)
							{
								adult[num_adults_chosen].maternal[jj].allelic_effects[qtls_locations[jj][jjj]] = progeny[j].maternal[jj].allelic_effects[qtls_locations[jj][jjj]];
								adult[num_adults_chosen].paternal[jj].allelic_effects[qtls_locations[jj][jjj]] = progeny[j].paternal[jj].allelic_effects[qtls_locations[jj][jjj]];
							}
						}
					}
					adult[num_adults_chosen].phenotype = progeny[j].phenotype;
					adult[num_adults_chosen].genotype = progeny[j].genotype;
					if (progeny[j].female){
						adult[num_adults_chosen].female = true;
						num_females++;
					}
					else{
						adult[num_adults_chosen].female = false;
						num_males++;
					}
					car_cap_rem = car_cap_rem - 1;
					num_adults_chosen++;
				}//end of if keep_prob
				else
					progeny[j].alive = false;
			}//end of if Alive
			prog_left = prog_left - 1;
		}//end of for j
		population_size = num_adults_chosen;
		if (population_size == 0)
			extinct = true;
	}//end Density Regulation

	void UpdateDemographics()
	{
		int j;
		num_males, num_females, mean_mal_trait = mean_fem_trait = 0;
		population_size = adult.size();
		for (j = 0; j < population_size; j++)
		{
			if (adult[j].female)
			{
				num_females++;
				mean_fem_trait = mean_fem_trait + adult[j].phenotype;
			}
			else
			{
				num_males++;
				mean_mal_trait = mean_mal_trait + adult[j].phenotype;
			}
		}
		mean_mal_trait = mean_mal_trait / num_males;
		mean_fem_trait = mean_fem_trait / num_females;
	}

	vector<bool> SamplePop(int sample_size, int pop_size)
	{
		//Pull a random sample of adults and a sample of offspring
		//Sort adults into males and females
		//Calculate allele frequencies for males, females, and offspring
		//Compare using Fst approaches
		if (sample_size == 0 || sample_size > pop_size)
			sample_size = pop_size;
		vector<bool> sampled_inds;
		int j, rand1, num_sampled;

		for (j = 0; j < pop_size; j++)
			sampled_inds.push_back(false);
		num_sampled = 0;
		while (num_sampled < sample_size)
		{
			rand1 = randnum(population_size);//need to sample without replacement
			if (sampled_inds[rand1] == false)
			{
				sampled_inds[rand1] = true;
				num_sampled++;
			}
		}
		return sampled_inds;
	}
};//population

void migration(population &pop_receiving, population &pop_giving, parameters &p)
{
	int j, jj, jjj, num_migrating;
	int rand1, rand2, rand3, rand4;

	num_migrating = pop_giving.adult.size() * p.migration_rate;

	while (num_migrating > 0)
	{
		rand1 = randnum(pop_giving.adult.size());
		pop_receiving.adult.push_back(pop_giving.adult[rand1]);
		pop_giving.adult.erase(pop_giving.adult.begin() + rand1);
		num_migrating--;
	}
	pop_receiving.UpdateDemographics();
	pop_giving.UpdateDemographics();
}

void output_genepop(vector<population> &pops, vector<vector<bool>> &sampled, parameters &p, string genepop_name)
{
	int j, jj, jjj, k;
	ofstream genepop;
	genepop.open(genepop_name);
	cout << "\ngenepop output file " << genepop_name << " open.\n";
	genepop << "Individual-based simulated data";
	for (j = 0; j < p.num_chrom; j++)
	{
		for (jj = 0; jj < p.num_markers; jj++)
		{
			genepop << "\nChrom" << j << "SNP" << jj;
		}
	}
	for (k = 0; k < p.num_pops; k++)
	{
		genepop << "\nPop";
		for (j = 0; j < pops[k].population_size; j++)
		{
			if (sampled[k][j])//check to see if the individual is sampled
			{
				genepop << "\nInd" << j << "\t,\t";
				for (jj = 0; jj < p.num_chrom; jj++)
				{
					for (jjj = 0; jjj < p.num_markers; jjj++)
					{
						//0101 indicates homozygous for allele 01
						if (p.num_alleles < 10)
						{
							genepop << "\t0" << pops[k].adult[j].maternal[jj].loci[jjj] << "0" << pops[k].adult[j].paternal[jj].loci[jjj];
						}
						if (p.num_alleles > 10)
						{
							genepop << '\t' << pops[k].adult[j].maternal[jj].loci[jjj] << pops[k].adult[j].paternal[jj].loci[jjj];
						}
					}
				}
			}
		}
	}//end of pops
	genepop.close();
}