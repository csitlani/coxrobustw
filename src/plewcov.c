/*
 * Copyright (C) 2006 Filip Borowicz
 *
 * Edits to original lin.c code to create this new plewcov.c code:
 * Copyright (C) 2018 Colleen Sitlani
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
 * MA  02111-1307  USA
 */

/*
 * @parma the exp(beta) vector
 * @parma the vector of times
 * @parma the vector of censoring indicators
 * @parma the covariate matrix z
 * @parma the vector of exp(prev_beta'z) (prev_beta - beta obatined through previous step)
 * @parma the number of individuals
 * @parma the number of covariates
 * @parma the sampling weights
 * @parma the result
 *
 * @return void
 */
#include "coxrobust.h"

void plewcov(double *exp_zbeta, double *time, int *status, double *covar,
         double *prev_exp_zbeta, int *n_row, int *n_col,
         double *sweight, double *result) {

    int i, j, k;
    double tmp;
    double sum1;
    double sum2;
    double *tmp1;
    double *tmp2;
    double **z;
    double **res;

    tmp1  = (double *)R_alloc(*n_row, sizeof(double));
    tmp2  = (double *)R_alloc(*n_row, sizeof(double));
    z     = dmatrix(covar, *n_row, *n_col);
    res   = dmatrix(result, *n_row, *n_col);


    for (j=0; j<*n_col; j++) {
        for (i=0; i<*n_row; i++) {
            
            if ( status[i] != 0 ) {

               sum1 = 0;
               sum2 = 0;
               for (k=i; k<*n_row; k++) {
                   tmp = exp_zbeta[k] * sweight[k];
                   sum1 += tmp;
                   sum2 += z[j][k]*tmp;
               }

               if ( sum1 == 0 ) {
                  sum1 = 1.0;
               }

               res[j][i] = ( z[j][i] - sum2/sum1 );
               tmp1[i]   = sum1;
	       tmp2[i]   = sum2; 

            } else {

               res[j][i] = 0;
               tmp1[i] = 1;
               tmp2[i] = 1;

            }

        }

        for (i=0; i<*n_row; i++) {
            for (k=0; k<=i; k++) {
                res[j][i] -= status[k] * exp_zbeta[i] * sweight[k] * 
			(z[j][i]-tmp2[k]/tmp1[k]) / tmp1[k];
            }
        }
    
    }

}
