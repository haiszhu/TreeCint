#include "cgto.h"
#include <stdlib.h>

void GTOval_sph(int ngrids, int *shls_slice, int *ao_loc,
                double *ao, double *coord, uint8_t *non0table,
                int *atm, int natm, int *bas, int nbas, double *env);

void GTOval_sph_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    GTOval_sph(*ngrids, shls_slice, ao_loc, ao, coord, non0tab_uint8, atm, *natm, bas, *nbas, env);

}

void GTOval_h2o_ccpvtz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    GTOval_sph(*ngrids, shls_slice, ao_loc, ao, coord, non0tab_uint8, atm, *natm, bas, *nbas, env);

}

void GTOval_nh3_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}

void GTOval_nh3_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}

void GTOval_h2o_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}

void GTOval_h2o_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}


void GTOval_uracil_dimer_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}

void GTOval_uracil_dimer_aug_ccpvdz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }

}

void GTOval_sph_generic_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[ (size_t) *nbas ];
    memset(non0tab_uint8, 1, (size_t) *nbas * sizeof(uint8_t));
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int i = 0; i < *ngrids; i++) {
      double *current_coord = &coord[i * 3];
      double *current_ao = &ao[i * nao];
      GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }
}

// evaluation of ith basis
void GTOival_sph_generic_mwrap(int *i, int *ngrids, int *shls_slice, int *ao_loc,
  double *aoi, double *coord, int *non0tab,
  int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[ (size_t) *nbas ];
    memset(non0tab_uint8, 1, (size_t) *nbas * sizeof(uint8_t));
    int nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
    #pragma omp parallel for
    for (int j = 0; j < *ngrids; j++) {
      double *current_coord = &coord[j * 3];
      // double *current_ao = &ao[j * nao];
      // GTOval_sph(1, shls_slice, ao_loc, current_ao, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
      // aoi[j] = current_ao[*i-1];
      double *current_aoi = &aoi[j];
      GTOival_sph(*i, 1, shls_slice, ao_loc, current_aoi, current_coord, non0tab_uint8, atm, *natm, bas, *nbas, env);
    }
}

