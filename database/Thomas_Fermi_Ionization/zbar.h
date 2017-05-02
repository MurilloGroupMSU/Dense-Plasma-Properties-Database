/* $Id:  $  */
#ifndef ZBAR_H
#define ZBAR_H
typedef struct  ion_state_st {  int  index; double T, Z, n, A,B,C, nf, dnf, vol, zBar;}  ION_STATE;
void zBarFunc(int nIonType,ION_STATE *ion);
void zBarFunc2(int nspec, double T, double *Z, double *n, double *zBar);
#endif
