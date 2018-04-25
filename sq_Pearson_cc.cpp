/* Implementation of 'sq_Pearson_cc.h' */
#include "sq_Pearson_cc.h"

#include <stdlib.h>
#include <math.h>


double sq_Pearson_cc::cal_r_squared(Init_Type & engine, int ** gibbs_genotypes, int snp_id)
{
	double varx = this->varianceX(gibbs_genotypes, snp_id, engine);
	double varz = this->varianceZ(snp_id, engine);
	double covxz = this->covariance(gibbs_genotypes, snp_id, engine);
	double r2 = pow(covxz, 2) / (varx * varz);
	if (r2 > 0.9999999999) r2 = 1.0;
	return r2;
}

double sq_Pearson_cc::cal_dosage(int * gibbs_genotypes, int indiv_id, int snp_id)
{
	return this->expectation(gibbs_genotypes, indiv_id, snp_id);
}

double sq_Pearson_cc::varianceX(int ** gibbs_genotypes, int snp_id, Init_Type & engine)
{
	int indiv_n = engine.individuals;
	double var = 0.0;
	double part1 = 0.0;
	double part2 = 0.0;
	double posterior_probs[3];
	double sum = 0.0;
	for (int i = 0; i != indiv_n; ++i) {
		sum = gibbs_genotypes[i][0] + gibbs_genotypes[i][1] + gibbs_genotypes[i][2];
		for (int j = 0; j != 3; ++j) {
			posterior_probs[j] = double(gibbs_genotypes[i][j]) / sum;
		}
		double section1 = posterior_probs[1] + 2 * posterior_probs[2];
		double section2 = posterior_probs[1] + 4 * posterior_probs[2];
		part1 += section2;
		part2 += section1;
	}
	part1 = part1 / indiv_n;
	part2 = pow(part2, 2) * (pow((1.0 / indiv_n), 2));
	var = part1 - part2;
	return var;
}

double sq_Pearson_cc::varianceZ(int snp_id, Init_Type & engine)
{
	int indiv_n = engine.individuals;
	double var = 0.0;
	double part1 = 0.0;
	double part2 = 0.0;
	for (int i = 0; i != indiv_n; ++i) {
		double genotype = engine.haplotypes[2 * i + 0].haplotype[snp_id];
		genotype += engine.haplotypes[2 * i + 1].haplotype[snp_id];
		part1 += pow(genotype, 2);
		part2 += genotype;
	}
	part1 = part1 / indiv_n;
	part2 = pow(part2, 2) * (pow((1.0 / indiv_n), 2));
	var = part1 - part2;
	return var;
}

double sq_Pearson_cc::covariance(int ** gibbs_genotypes, int snp_id, Init_Type & engine)
{
	int indiv_n = engine.individuals;
	double cov = 0.0;
	double part1 = 0.0;
	double part2 = 0.0;
	double part21 = 0.0;
	double part22 = 0.0;
	double posterior_probs[3];
	double sum = 0.0;
	for (int i = 0; i != indiv_n; ++i) {
		sum = gibbs_genotypes[i][0] + gibbs_genotypes[i][1] + gibbs_genotypes[i][2];
		for (int j = 0; j != 3; ++j)
			posterior_probs[j] = double(gibbs_genotypes[i][j]) / sum;
		double genotype = engine.haplotypes[2 * i + 0].haplotype[snp_id];
		genotype += engine.haplotypes[2 * i + 1].haplotype[snp_id];
		double section = posterior_probs[1] + 2 * posterior_probs[2];
		part1 += (genotype * section);
		part21 += section;
		part22 += genotype;
	}
	part1 = part1 / indiv_n;
	part2 = (pow((1.0 / indiv_n), 2)) * part21 * part22;
	cov = part1 - part2;
	return cov;
}

double sq_Pearson_cc::expectation(int * gibbs_genotypes, int indiv_id, int snp_id)
{
	double dosage = 0.0;
	double posterior_probs[3];
	double sum = gibbs_genotypes[0] + gibbs_genotypes[1] + gibbs_genotypes[2];
	for (int g = 0; g <= 2; ++g) {
		posterior_probs[g] = double(gibbs_genotypes[g]) / sum;
		dosage += posterior_probs[g] * g;
	}
	return dosage;
}


/* END */
