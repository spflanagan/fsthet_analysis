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
	int end, i, ii, iii,r, rand_pop, ld_count;
	double dp, avg_dp;
	vector<population> pops;
	vector<int> qtls;
	parameters params;
	params.set_parameters();
	
	string ld_file_name = params.base_file_name + "ld.txt";
	ofstream ld_file;
	ld_file.open(ld_file_name);

	for (r = 0; r < params.num_reps; r++)
	{
		ld_file << "Rep " << r;
		for (ii = 0; ii < (params.num_chrom*params.num_qtls); ii++)
		{
			qtls.push_back(randnum(params.num_chrom*params.num_markers));
		}
		sort(qtls.begin(), qtls.end());

		for (i = 0; i < params.num_pops; i++)
		{//establish the populations
			pops.push_back(population());
			cout << "\nInitializing Pop " << i;
			pops[i].initialize(qtls, params);
		}

		for (i = 0; i < params.generations; i++)
		{//each generation there is a chosen one. she alone will stand against the vampires, the demons, and the forces of darkness. she is the slayer.
			ld_file << '\n' << i;
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
				//calculate ld
				avg_dp = 0;
				for (int c = 0; c < params.num_chrom; c++)
				{
					ld_count = 0;
					while (ld_count < 100)
					{//random within-chromosomes
						dp = pops[ii].calc_ld(c, c, randnum(params.num_markers), randnum(params.num_markers), params);
						if (dp != -5)
						{
							avg_dp = avg_dp + dp;
							ld_count++;
						}
					}
				}
				avg_dp = avg_dp / (100 * params.num_chrom);
				ld_file << '\t' << avg_dp;
			}//end numpops
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