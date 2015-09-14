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

using namespace std;

class parameters
{
public:
	int pop_size, carrying_capacity;
	int max_fecundity, max_encounters, max_progeny;
	int total_num_loci, num_poly_loci;
	int sample_size, generations, num_pops;
	int num_chrom, num_markers, num_qtls, num_alleles;
	double recomb_rate, mut_rate, mut_var, env_var;
	double GaussianPrefMean, migration_rate;

	void set_parameters()
	{
		//Initialize Parameters
		carrying_capacity = 5000;
		max_fecundity = 4;
		max_progeny = carrying_capacity*max_fecundity;
		num_chrom = 4;
		num_markers = 1000;
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
		sample_size = 4000;//set sample size here
		generations = 2000;
		num_pops = 10;
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
	vector<vector<bool>> qtls_locations;
	int population_size, num_males, num_females, num_progeny;
	double mean_mal_trait, mean_fem_trait;

public:
	void initialize(vector <int> &qtls, parameters &p)
	{
		int j,jj,jjj, count,q;
		//set up qtl tracker
		count = q = 0;
		for (j = 0; j < p.num_chrom; j++)
		{
			qtls_locations.push_back(vector<bool>());
			for (jj = 0; jj < p.num_markers; jj++)
			{
				if (count == qtls[q])
				{
					qtls_locations[j].push_back(true);
					q++;
				}
				else
					qtls_locations[j].push_back(false);
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
				for (jjj = 0; jjj < p.num_qtls; jjj++)
				{
					adult[j].maternal[jj].allelic_effects[jjj] = temp_allele[j%p.num_alleles];
					adult[j].paternal[jj].allelic_effects[jjj] = temp_allele[j%p.num_alleles];
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
						if (qtls_locations[WhichChromosome][RCj])
							RecombinedChr.allelic_effects[RCj] = Parent.maternal[WhichChromosome].allelic_effects[RCj];
					}
				}
				else
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
					{
						RecombinedChr.loci[RCj] = Parent.paternal[WhichChromosome].loci[RCj];
						if (qtls_locations[WhichChromosome][RCj])
							RecombinedChr.allelic_effects[RCj] = Parent.paternal[WhichChromosome].allelic_effects[RCj];
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
				for (RCi = 0; RCi < p.num_qtls; RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.maternal[WhichChromosome].allelic_effects[RCi];
			}
			else
			{
				for (RCi = 0; RCi < p.num_markers; RCi++)
					RecombinedChr.loci[RCi] = Parent.paternal[WhichChromosome].loci[RCi];
				for (RCi = 0; RCi < p.num_qtls; RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.paternal[WhichChromosome].allelic_effects[RCi];
			}
		}//else (no recomb)
	}//end RecombineChromosome

	void Mating(parameters &p)
	{
		int cc;
		int NumProg;
		double dubrand;
		NumProg = 0;
		//Mate Choice:
		int mm, nn, males;
		int Encounters, MaleID, irndnum, Females;
		double MeanFemaleTrait;
		bool MateFound, MalesPresent;
		double MeanMaleTrait, MateProb;
		int MaleIndex, FemaleID;
		int Counter3;
		vector <int> MaleList;
		NumProg = 0;
		Mated = 0;
		Females = 0;
		ProgFem = 0;
		ProgMale = 0;
		NumFemales = 0;
		MalesPresent = false;
		NumMales = 0;
		MeanMaleTrait = 0;
		double SDMaleTrait = 0;
		MeanFemaleTrait = 0;
		Counter1 = 0;
		Counter2 = 0;
		Counter3 = 0;
		//determine mean male trait
		for (males = 0; males < PopulationSize; males++)
		{
			Adult[males].MateFound = 0;
			if (!Adult[males].IsFemale)
			{
				MalesPresent = true;
				MaleList.push_back(males);
				NumMales++;
				MeanMaleTrait = MeanMaleTrait + Adult[males].phenotype;
			}
		} // end of males
		if (NumMales > 0) {
			MeanMaleTrait = MeanMaleTrait / NumMales;
		}
		else {
			MeanMaleTrait = 0;
			popExtinct = true;
		}
		for (males = 0; males < PopulationSize; males++){
			if (!Adult[males].IsFemale)
				SDMaleTrait = SDMaleTrait + (Adult[males].phenotype - MeanMaleTrait)*(Adult[males].phenotype - MeanMaleTrait);
		}
		SDMaleTrait = sqrt(SDMaleTrait / NumMales);
		for (mm = 0; mm < PopulationSize; mm++)
		{
			MateFound = false;
			if (Adult[mm].IsFemale && MalesPresent)
			{
				FemaleID = mm;
				MeanFemaleTrait = MeanFemaleTrait + Adult[mm].phenotype;
				Females++;
				Encounters = 0;
				while (!MateFound && Encounters <= MaximumEncounters)
				{
					irndnum = randnum(NumMales);
					MaleIndex = MaleList[irndnum];
					if (GaussianPrefVariance > 0)
						MateProb = exp(-0.5 * (Adult[MaleIndex].phenotype - GaussianPrefMean)*
						(Adult[MaleIndex].phenotype - GaussianPrefMean) / GaussianPrefVariance);
					else
						MateProb = 1;
					dubrand = genrand();
					if (dubrand < MateProb)
					{
						MateFound = true;
						MaleID = MaleIndex;
						Fathers[Counter2] = MaleID;
						Adult[MaleID].MateFound++;
						Counter2++;
						Mated++;
					}
					Encounters++;
				}//while
				if (MateFound)
				{
					Mothers[Counter1] = mm;
					Counter1++;
					Adult[mm].MateFound++;
					if (NumProg >= MaxNumProg)
						NumProg = MaxNumProg - 1;
					//mother is Parent 1, mm
					//father is parent 2, MateID
					for (nn = 0; nn < MaximumFecundity; nn++)
					{
						Progeny[NumProg].Alive = true;
						//calculate phenotype in mutation once the genotype is determined
						for (cc = 0; cc < NumChrom; cc++)//go through chromosome by chromosome
						{
							RecombineChromosome(Progeny[NumProg].maternal[cc], Adult[FemaleID], cc, RecombRate);
							RecombineChromosome(Progeny[NumProg].paternal[cc], Adult[MaleID], cc, RecombRate);
						}//end of chromosome
						if (genrand() < 0.5)
						{
							Progeny[NumProg].IsFemale = true;
							ProgFem++;
						}
						else
						{
							Progeny[NumProg].IsFemale = false;
							ProgMale++;
						}
						NumProg++;
					}//for nn
				}
				if (!MateFound)
				{//Keep track of the females that didn't mate
					Nonmated[Counter3] = mm;
					Counter3++;
				}
			}//if female
		}//end of mm
		ProgenyNum = NumProg;
		NumFemales = Females;
		NumNotMated = Counter3;
		MeanFemaleTrait = MeanFemaleTrait / Females;
	}//end mating

	void Mutation()
	{
		int m, mm, irand, irand2, gg, ggg, irand3, locus;
		double rnd1, rnd2;
		double IndMutationRate;
		double MutSD;
		bool mutated;

		MutSD = sqrt(MutationalVariance);
		IndMutationRate = MutationRate * 2 * TotalLociNo;

		for (m = 0; m < ProgenyNum; m++)
		{
			rnd1 = genrand();
			Progeny[m].phenotype = 0;
			Progeny[m].Genotype = 0;
			mutated = false;
			if (rnd1 < IndMutationRate)
			{
				irand = randnum(NumChrom);
				irand2 = randnum(NumMarkers);
				locus = MarkerLoci[irand].LociOnChrom[irand2];
				rnd2 = genrand();//to choose maternal or paternal
				if (rnd2 < 0.5)//affects maternal chromosome	
				{
					while (!mutated){
						irand3 = randnum(NumAlleles);
						if (!Progeny[m].maternal[irand].LociArray[locus] == irand3)
						{
							Progeny[m].maternal[irand].LociArray[locus] = irand3;
							mutated = true;
						}
					}
					for (mm = 0; mm < NumQTLs; mm++)
					{
						if (Locations[irand].LociOnChrom[mm] == irand2)
							Progeny[m].maternal[irand].allelicEffects[mm] =
							Progeny[m].maternal[irand].allelicEffects[mm] + randnorm(0, MutSD);
					}
				}
				else//affects paternal chromosome
				{
					while (!mutated){
						irand3 = randnum(NumAlleles);
						if (!Progeny[m].paternal[irand].LociArray[locus] == irand3)
						{
							Progeny[m].paternal[irand].LociArray[locus] = irand3;
							mutated = true;
						}
					}
					for (mm = 0; mm < NumQTLs; mm++)
					{
						if (Locations[irand].LociOnChrom[mm] == irand2)
							Progeny[m].paternal[irand].allelicEffects[mm] =
							Progeny[m].paternal[irand].allelicEffects[mm] + randnorm(0, MutSD);
					}
				}
			}//end of if

			for (gg = 0; gg < NumChrom; gg++)
			{
				for (ggg = 0; ggg < NumQTLs; ggg++){
					Progeny[m].Genotype = Progeny[m].Genotype +
						Progeny[m].maternal[gg].allelicEffects[ggg] + Progeny[m].paternal[gg].allelicEffects[ggg];
				}
			}
			Progeny[m].phenotype = Progeny[m].Genotype + randnorm(0, EnvStdDev);
		}//end of m
	}//mutation

	void SelectionOnPhenotypes(double dSelectionStrength)
	{
		int i, ProgAlive;
		double dSurvProb;
		double drnum1;
		double dOptimum;
		double phenSD = 0;
		double phenMean = 0;
		double num;
		int malecount = 0;
		ProgAlive = 0;
		//calc mean male phenotype & std dev
		for (i = 0; i < ProgenyNum; i++){
			if (!Progeny[i].IsFemale){
				phenMean = phenMean + Progeny[i].phenotype;
				malecount++;
			}
		}
		num = malecount;
		phenMean = phenMean / num;
		for (i = 0; i < ProgenyNum; i++){
			if (!Progeny[i].IsFemale)
				phenSD = phenSD + (Progeny[i].phenotype - phenMean)*(Progeny[i].phenotype - phenMean);
		}
		phenSD = sqrt(phenSD / num);

		dOptimum = 0;
		for (i = 0; i < ProgenyNum; i++) {
			if (Progeny[i].IsFemale)
			{
				Progeny[i].Alive = true;
				ProgAlive++;
			}
			else//selection only on males
			{
				if (dSelectionStrength > 0)
					dSurvProb = exp(-1 * (Progeny[i].phenotype - dOptimum)*(Progeny[i].phenotype - dOptimum)
					/ (2 * dSelectionStrength));
				else
					dSurvProb = 1;
				//cout<<dSurvProb<<'\n';
				drnum1 = genrand();
				if (drnum1 < dSurvProb)
				{
					Progeny[i].Alive = true;
					ProgAlive++;
				}
				else
					Progeny[i].Alive = false;
			}
		} // end of i
	}//end selection on phenotypes

	void DensityRegulation()
	{
		int p, pp, ppp;
		int iNumAdultsChosen;
		double CarCapUnfilled, ProgLeft, KeepProb;
		double DRrandnum;
		NumMales = 0;
		NumFemales = 0;
		ProgLeft = 0;
		//count the ones that are still alive
		for (p = 0; p < ProgenyNum; p++)
		{
			if (Progeny[p].Alive)
				ProgLeft++;
		}
		CarCapUnfilled = CarryingCapacity;
		iNumAdultsChosen = 0;
		for (p = 0; p<ProgenyNum; p++)
		{
			if (Progeny[p].Alive)
			{
				if (ProgLeft == 0)
					KeepProb = 0;
				else
					KeepProb = CarCapUnfilled / ProgLeft;
				DRrandnum = genrand();
				if (DRrandnum<KeepProb)
				{//then turn it into an adult
					Adult[iNumAdultsChosen].Alive = true;
					Adult[iNumAdultsChosen].MateFound = 0;
					for (pp = 0; pp < NumChrom; pp++)
					{
						for (ppp = 0; ppp < NumMarkers; ppp++)
						{
							Adult[iNumAdultsChosen].maternal[pp].LociArray[ppp] = Progeny[p].maternal[pp].LociArray[ppp];
							Adult[iNumAdultsChosen].paternal[pp].LociArray[ppp] = Progeny[p].paternal[pp].LociArray[ppp];
						}
						for (ppp = 0; ppp < NumQTLs; ppp++)
						{
							Adult[iNumAdultsChosen].maternal[pp].allelicEffects[ppp] = Progeny[p].maternal[pp].allelicEffects[ppp];
							Adult[iNumAdultsChosen].paternal[pp].allelicEffects[ppp] = Progeny[p].paternal[pp].allelicEffects[ppp];
						}
					}
					Adult[iNumAdultsChosen].phenotype = Progeny[p].phenotype;
					Adult[iNumAdultsChosen].Genotype = Progeny[p].Genotype;
					if (Progeny[p].IsFemale){
						Adult[iNumAdultsChosen].IsFemale = true;
						ChosenSexes[iNumAdultsChosen] = 1;
						NumFemales++;
					}
					else{
						Adult[iNumAdultsChosen].IsFemale = false;
						ChosenSexes[iNumAdultsChosen] = 0;
						NumMales++;
					}
					CarCapUnfilled = CarCapUnfilled - 1;
					iNumAdultsChosen++;
				}//end of if KeepProb
				else
					Progeny[p].Alive = false;
			}//end of if Alive
			ProgLeft = ProgLeft - 1;
		}//end of for p
		PopulationSize = iNumAdultsChosen;
		if (PopulationSize == 0)
			popExtinct = true;
	}//end Density Regulation

};//population
