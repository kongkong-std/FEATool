#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <cjson/cJSON.h>
#include <stdbool.h>
#include <parmetis.h>
#include <metis.h>

// function prototype
int LinsysDiagScal(Mat *mat /*diagonal scaling matrix*/, Vec *rhs /*diagonal scaling r.h.s.*/, Vec *diag /*processed diagonal*/);
int SolDiagScal(Vec *tilde_x /*tilde solution*/, Vec *diag /*processed diagonal*/, Vec *x /*solution*/);

#endif // main.h