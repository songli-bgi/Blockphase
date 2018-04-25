#include "EM_hfs.h"

#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#define MAX_ITER 100

EM_hfs::~EM_hfs()
{
	int size = haps_.size();
	for (int i = 0; i != size; ++i) {
		delete haps_[i];
	}
}

bool EM_hfs::clean_array(uint8_t * pArr, int len, int kind)
{
	for (int i = len - 1; i >= 0; --i) {
		if (i == 0) {
			pArr[i] = 1;
			break;
		}
		if (pArr[i] / kind == 1) {
			pArr[i - 1] += 1;
			pArr[i] %= kind;
		} else {
			break;
		}
	}
	return pArr[0] > 0 ? false : true;
}

int EM_hfs::enumerate(int snp_n)
{
	snp_n_ = snp_n;
	int idx = 0;
	uint8_t * pArr = new uint8_t [snp_n_ + 1];
	for (int i = 0; i != snp_n_ + 1; ++i)
		pArr[i] = 0;
	bool bLoop = true;
	int kind = 2; // 0 or 1
	while (bLoop) {
		bLoop = this->clean_array(pArr, snp_n_ + 1, kind);
		if (bLoop) {
			Haplotype * hap_tmp = new Haplotype(snp_n_);
			memcpy(hap_tmp->haplotype, pArr + 1, snp_n_ * sizeof(*pArr));
			haps_.push_back(hap_tmp);
			freqs_[*hap_tmp] = 0.;
			hap_idx_[*hap_tmp] = idx++;
			pArr[snp_n_]++;
		}
	}
	delete [] pArr;
	return haps_.size();
}


void EM_hfs::initialize(_GenoHap & genehap, hap_t * haplotypes, int first_p, int site_n)
{
	for (int beg = 0; beg < site_n - 2; ) {
		if (genehap.whole_reference == NULL) this->EM_iterate(genehap, first_p, beg);
		else { // initialize using reference panel
			this->Ref_estimate(genehap, first_p, beg);
		}
		this->init_haps(genehap, first_p, beg, haplotypes);
		if (beg + snp_n_ < site_n - snp_n_) beg += (snp_n_ - 1);
		else if (beg == site_n - snp_n_) break;
		else beg = site_n - snp_n_;
	}
	hap2geno(genehap);
	return;
}

void EM_hfs::hap2geno(_GenoHap & genohap)
{
	int individuals = genohap.genotypes.size();
	int snp_number = genohap.genotypes[ 0 ].size();
	for ( int i = 0; i < individuals; i++ ) {
		for ( int j = 0; j < snp_number; j++ ) {
			genohap.genotypes[ i ][ j ] = genohap.haplotypes[ i ].first[ j ] + genohap.haplotypes[ i ].second[ j ];
		}
	}
}

double EM_hfs::cal_likelihood(const _GenoHap & genehap, int idv, int first_p, int beg, int h1_idx, int h2_idx)
{
	double lh = 1.;
	for (int i = 0; i != snp_n_; ++i) {
		uint8_t g = haps_[h1_idx]->haplotype[i] + haps_[h2_idx]->haplotype[i];
		lh *= genehap.glf[idv * genehap.snp_number * 3 + (first_p + beg + i) * 3 + g];
	}
	return lh;
}

double EM_hfs::cal_prior(int h1_idx, int h2_idx)
{
	double prior = 0.;
	if (h1_idx == h2_idx) {
		prior = freqs_[*haps_[h1_idx]] * freqs_[*haps_[h2_idx]];
	} else {
		prior = 2 * freqs_[*haps_[h1_idx]] * freqs_[*haps_[h2_idx]];
	}
	return prior;
}

void EM_hfs::EM_iterate(const _GenoHap & genehap, int first_p, int beg)
{
	int size = haps_.size();
	vector<double> N_c(size, 0.);
	double p = 1. / (size + 0.);
	for (map<Haplotype, double>::iterator it = freqs_.begin(); it != freqs_.end(); ++it) {
		it->second = p;
	}
	for (int iter = 0; iter < MAX_ITER; ++iter) {
		for (int i = 0; i != genehap.individuals; ++i) {
			vector<double> N_cc(size, 0.);
			double t_cc = 0.;
			for (int h1 = 0; h1 < size; ++h1) {
				for (int h2 = h1; h2 < size; ++h2) {
					double lh = this->cal_likelihood(genehap, i, first_p, beg, h1, h2);
					lh *= this->cal_prior(h1, h2);
					N_cc[h1] += lh / 2;
					N_cc[h2] += lh / 2;
					t_cc += lh;
				}
			}
			for (int h = 0; h < size; ++h) {
				N_c[h] += (N_cc[h] / t_cc);
			}
		} // E step;
		double tot = 0.;
		for (vector<double>::iterator it = N_c.begin(); it != N_c.end(); ++it) {
			tot += *it;
		}
		for (int i = 0; i != size; ++i) {
			N_c[i] /= tot;
		} // M step
		double norm = 0.;
		for (int i = 0; i != size; ++i) {
			norm += pow((N_c[i] - freqs_[*haps_[i]]), 2);
			freqs_[*haps_[i]] = N_c[i];
		}
		if (norm < 1e-6) break;
	}
	return;
}


void EM_hfs::Ref_estimate(const _GenoHap & genehap, int first_p, int beg)
{
	int size = haps_.size();
	vector<int> h_c(size, 0);
	Haplotype hap(snp_n_);
	for (int i = 0; i != genehap.whole_n_reference; ++i) {
		memcpy(hap.haplotype, genehap.whole_reference[i].haplotype + first_p + beg, snp_n_ * sizeof(*hap.haplotype));
		h_c[hap_idx_[hap]]++;
	}
	for (int i = 0; i != size; ++i) {
		double p = (h_c[i] + 0.) / (genehap.whole_n_reference * 2.);
		freqs_[*haps_[i]] = p;
	}
	return;
}


void EM_hfs::init_haps(const _GenoHap & genehap, int first_p, int beg, hap_t * haplotypes)
{
	double max_p = 0.;
	int size = haps_.size();
	int h1_idx, h2_idx;
	for (int i = 0; i != genehap.individuals; ++i) {
		double max_p = 0.;
		int h1_idx, h2_idx;
		h1_idx = h2_idx = -1;
		for (int h1 = 0; h1 < size; ++h1) {
			for (int h2 = h1; h2 < size; ++h2) {
				double lh = this->cal_likelihood(genehap, i, first_p, beg, h1, h2);
				double prior = this->cal_prior(h1, h2);
				double post = lh * prior;
				if (post > max_p) {
					max_p = post;
					h1_idx = h1;
					h2_idx = h2;
				}
			}
		}
		if (beg == 0 || haplotypes[i].first[beg] == haplotypes[i].second[beg]) {
			for (int j = 0; j != snp_n_; ++j) {
				haplotypes[i].first[beg + j] = haps_[h1_idx]->haplotype[j];
				haplotypes[i].second[beg + j] = haps_[h2_idx]->haplotype[j];
			}
		} else {
			if (haplotypes[i].first[beg] == haps_[h2_idx]->haplotype[0] && haplotypes[i].second[beg] == haps_[h1_idx]->haplotype[0]) {
				for (int j = 1; j != snp_n_; ++j) {
					haplotypes[i].first[beg + j] = haps_[h2_idx]->haplotype[j];
					haplotypes[i].second[beg + j] = haps_[h1_idx]->haplotype[j];
				}
			} else {
				for (int j = 1; j != snp_n_; ++j) {
					haplotypes[i].first[beg + j] = haps_[h1_idx]->haplotype[j];
					haplotypes[i].second[beg + j] = haps_[h2_idx]->haplotype[j];
				}
			}
		}
	}
	return;
}


void EM_hfs::display()
{
	int size = haps_.size();
	for (int i = 0; i != size; ++i) {
		cerr << i << " : ";
		for (int j = 0; j != snp_n_; ++j) {
			cerr << int(haps_[i]->haplotype[j]);
		}
		cerr << endl;
	}
}


