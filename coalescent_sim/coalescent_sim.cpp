//Author: Sarah P. Flanagan
//Date: 11 September 2015
//Purpose: write a coalescent simulation like FDIST and Lositan use to compare individual-based and coalescent methods
//Do coalescent simulations over-simplify to the extent that the null model doesn't do a good job of approximating the actual distribution?

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include "coalescent_classes.h"
#include "random_numbers.h"

using namespace std;

int main()
{
	struct node *tree;
	int sample_size, number_samples, in, *spectrum, nmuts, ndes, node;
	double theta, time, total_muts, time_of_size_change, factor;

	scanf_s(" %lf", &theta);
	scanf_s(" %d", &sample_size);
	scanf_s(" %d", &number_samples);
	scanf_s(" %lf", &time_of_size_change);
	scanf_s(" %lf", &factor);

	tree = (struct node *)malloc(2 * sample_size*sizeof(struct node));
	spectrum = (int *)malloc(sample_size*sizeof(int));

	for (in = 0; in < sample_size; in++) spectrum[in] = 0;
	total_muts = 0.0;

	for (in = 0; in < number_samples; in++)
	{
		make_tree(tree, sample_size);
		bottleneck(tree, sample_size, time_of_size_change, factor);
		for (node = 0; node < sample_size * 2 - 2; node++)
		{
			time = (tree[node].ancestor->time) - tree[node].time;
			nmuts = poissonrand(time*theta / 2.0);

			if (nmuts > 0)
			{
				ndes = count_desc(tree + node);
				spectrum[ndes] += nmuts;
				total_muts += nmuts;
			}
		}
	}

	printf(" Average number of mut's per sample: %lf\n\n", total_muts / number_samples);
	printf("freq of mut\tfraction of mutations\n\n");
	for (in = 1; in < sample_size; in++)
		printf("%d\t%lf\n", in, spectrum[in] / total_muts);
	cout << "input integer to quit:\n";
	cin >> in;
	return in;
}