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
	int end, i, ii, iii;
	vector<population> pops;
	vector<int> qtls;
	parameters params;
	params.set_parameters();


	for (ii = 0; ii < (params.num_chrom*params.num_qtls); ii++)
	{
		qtls.push_back(randnum(params.num_chrom*params.num_markers));
	}
	for (i = 0; i < params.num_pops; i++)
	{//establish the populations
		pops.push_back(population());
		pops[i].initialize(qtls, params);
	}

	for (i = 0; i < params.generations; i++)
	{//each generation there is a chosen one. she alone will stand against the vampires, the demons, and the forces of darkness. she is the slayer.
		//mating (includes recombination)

		//mutation

	}
	

	cout << "\ninput integer to quit\n";
	cin >> end;
	return 0;
}