#ifndef CGTO_H
#define CGTO_H

#include <stdint.h>
#include <string.h>

void GTOval_sph(int ngrids, int *shls_slice, int *ao_loc,
                double *ao, double *coord, uint8_t *non0table,
                int *atm, int natm, int *bas, int nbas, double *env);

void GTOival_sph(int i, int ngrids, int *shls_slice, int *ao_loc,
                double *aoi, double *coord, uint8_t *non0table,
                int *atm, int natm, int *bas, int nbas, double *env);

void GTOval_sph_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_h2o_ccpvtz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_nh3_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_nh3_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_h2o_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_h2o_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_uracil_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_uracil_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_sph_generic_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOival_sph_generic_mwrap(int *i, int *ngrids, int *shls_slice, int *ao_loc,
                      double *aoi, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);                      

#endif // CGTO_H
