/**
 * @file 		 sq_Pearson_cc.h
 * @brief 		 calculate squared Pearson correlation coefficient
 * @data 		 25-Jul-2013
 * @author 		 Jin Wei(jinwei@techsolutions.com)
 *
 * Calculate imputation R-squared of each site
 *
 */

#ifndef __SQ_PEARSON_CC_H__
#define __SQ_PEARSON_CC_H__

#include <stdint.h>
#include "Init_Type.h"

/// squared Pearson correlation coefficient & R-squared
class sq_Pearson_cc
{
 public:
	sq_Pearson_cc() {};
	~sq_Pearson_cc() {};

	/// Calculate R-squared and dosage one by one
	double cal_r_squared(Init_Type & engine, int ** gibbs_genotypes, int snp_id);
	double cal_dosage(int * gibbs_genotypes, int indiv_id, int snp_id);

 private:
	/// Prevent copy and assignment
	sq_Pearson_cc(const sq_Pearson_cc & );
	sq_Pearson_cc & operator = (const sq_Pearson_cc & );

	/// Variance of genotype posterior probability
	double varianceX(int ** gibbs_genotypes, int snp_id, Init_Type & engine);

	/// Variance of final decided genotypes
	double varianceZ(int snp_id, Init_Type & engine);

	/// Covariance of genotypes and genotype probability
	double covariance(int ** gibbs_genotypes, int snp_id, Init_Type & engine);

	/// Expectation of genotypes or dosage
	double expectation(int *  gibbs_genotypes, int indiv_id, int snp_id);

};


#endif /* sq_Pearson_cc.h */
