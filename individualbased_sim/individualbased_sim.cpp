#include "random_numbers.h"
#include "population.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <array>

using namespace std;

int main()
{
	int end, i, ii, iii,r, rand_pop;
	vector<population> pops;
	vector<int> qtls;
	parameters params;
	params.set_parameters();

	for (r = 0; r < params.num_reps; r++)
	{
		for (ii = 0; ii < (params.num_chrom*params.num_qtls); ii++)
		{
			qtls.push_back(randnum(params.num_chrom*params.num_markers));
		}
		for (i = 0; i < params.num_pops; i++)
		{//establish the populations
			pops.push_back(population());
			cout << "\nInitializing Pop " << i;
			pops[i].initialize(qtls, params);
		}

		for (i = 0; i < params.generations; i++)
		{//each generation there is a chosen one. she alone will stand against the vampires, the demons, and the forces of darkness. she is the slayer.
			for (ii = 0; ii < params.num_pops; ii++)
			{
				//mating (includes recombination)
				pops[ii].Mating(params);
				//mutation
				pops[ii].Mutation(params);
				//selection
				pops[ii].SelectionOnPhenotypes(params.viability_sel);
				//survival to adulthood
				pops[ii].DensityRegulation(params);
				//migration
				rand_pop = ii;
				while (rand_pop == ii)
					rand_pop = randnum(params.num_pops);
				migration(pops[ii], pops[rand_pop], params);
			}
			if ((i + 1) % 100 == 0)
				cout << "\nGeneration " << i + 1 << " complete.\n";
		}

		vector<vector<bool>> sampled_inds;
		for (i = 0; i < params.num_pops; i++)
		{
			sampled_inds.push_back(pops[i].SamplePop(params.sample_size, pops[i].population_size));
		}
		cout << "\nWriting individual genotypes to genepop input file.\n";
		stringstream genepop_out_name;
		genepop_out_name << "sim_genepop_rep" << r << ".txt";
		output_genepop(pops, sampled_inds, params, genepop_out_name.str());

		cout << "\nRep " << r << " complete.\n";
	}
	cout << "\ninput integer to quit\n";
	cin >> end;
	return 0;
}