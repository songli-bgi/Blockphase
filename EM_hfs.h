/**
 * @brief EM algorithm to estimate haplotype frequency spectrum
 * @author Jinwei (wei.jin@bgitechsolutions.com)
 * @data 1-14-2014
 * @version 1.0
 */
#ifndef __EM_HFS_H__
#define __EM_HFS_H__

#include <map>

#include "GenoHap_Init.h"
#include "Haplotype.h"

using std::map;
using std::pair;

struct Haplotype; // haplotype data structor
struct _GenoHap;  // store individual haplotypes and reference panel
struct hap_t;

struct hapcomp 
{
	bool operator () (const Haplotype & lhs, const Haplotype & rhs)
	{
		if (lhs.length < rhs.length) {
			return 1;
		} else if (lhs.length > rhs.length) {
			return 0;
		} else {
			if (memcmp(lhs.haplotype, rhs.haplotype, lhs.length) < 0) {
				return 1;
			} else {
				return 0;
			}
		}
		return 0;
	}
};

class EM_hfs
{
 public:
	/// Constructor
	EM_hfs() {};
	/// Destructor
	~EM_hfs();

	/// Enumerate all possible haplotype under the condition of given bi-allelic sites' number
	/**
	 * @param snp_n 	 snp number
	 * @return 			 number of possible haplotypes
	 */
	int enumerate(int snp_n);

	/// Haplotype initialization
	/**
	 * @param genehap 	 struct in which GL information is
	 * @param haplotypes struct to store the haplotypes
	 * @param first_p 	 position of the first site
	 * @param site_n 	 number of sites in genehap
	 */
	void initialize(_GenoHap & genehap, hap_t * haplotypes, int first_p, int site_n);

	void hap2geno(_GenoHap & genehap);
	/// EM algorithm to estimate haplotype frequency
	/**
	 * @param genehap 	 struct in which GL information is
	 * @param first_p 	 position of the first site
	 * @param beg 		 begin position on the chromosome region
	 */
	
	void EM_iterate(const _GenoHap & genehap, int first_p, int beg);

	/// Use reference panel to estimate haplotype frequency
	/**
	 * @param genehap 	 struct in which GL information is
	 * @param first_p 	 position of the first site
	 * @param beg 		 begin position on the chromosome region
	 */
	void Ref_estimate(const _GenoHap & genehap, int first_p, int beg);

	/// Initialize haplotypes
	/**
	 * @param genehap 	 struct in which GL information is
	 * @param first_p 	 position of the first site
	 * @param beg 		 begin position on the chromosome region
	 * @param haplotypes initialized haplotypes
	 */
	void init_haps(const _GenoHap & genehap, int first_p, int beg, hap_t * haplotypes);

	/* Chech the result */
	/// Display all the haplotypes
	void display();


 private:
	/// Prevent copy and assignment
	EM_hfs(const EM_hfs & );
	EM_hfs & operator = (const EM_hfs & );

	/// Called by enumerate
	bool clean_array(uint8_t * pArr, int len, int kind);

	/// Calculated likelihood of given diplotype at specified individual
	double cal_likelihood(const _GenoHap & genehap, int idv, int first_p, int beg, int h1_idx, int h2_idx);

	/// Calculated the prior of haplotype pair
	double cal_prior(int h1_idx, int h2_idx);

	map<Haplotype, double, hapcomp> freqs_; 	 ///< haplotype frequency
	map<Haplotype, int, hapcomp> hap_idx_; 		 ///< haplotype index
	vector<Haplotype *> haps_; 					 ///< all possible haplotypes
	int snp_n_;

};


#endif /* EM_hfs.h */
