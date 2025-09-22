
#include <petscksp.h>

int main(int argc, char **argv)
{
  // PetscErrorCode ierr;
  Mat A = NULL;
  // KSP ksp;
  PetscViewer viewer;
  char fileA[PETSC_MAX_PATH_LEN] = "A.petsc";
  char output[PETSC_MAX_PATH_LEN];
  PetscBool flg;

  PetscInitialize(&argc, &argv, (char *)0, NULL);

  /* --- get filenames from command line options (if provided) --- */
  PetscOptionsGetString(NULL, NULL, "-A", fileA, sizeof(fileA), &flg);
  // PetscOptionsGetString(NULL,NULL,"-B",fileB,sizeof(fileB),&flg);
  PetscOptionsGetString(NULL, NULL, "-output", output, PETSC_MAX_PATH_LEN, &flg);

  /* --- load matrix A --- */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileA, FILE_MODE_READ, &viewer);
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetFromOptions(A);
  MatLoad(A, viewer);
  PetscViewerDestroy(&viewer);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  int nrows, ncols;
  MatGetSize(A, &nrows, &ncols);
  PetscPrintf(PETSC_COMM_WORLD, "-------row: %d\n", nrows);
  // MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==== basic information of linear system ====\n"));
  PetscInt m_mat = 0, n_mat = 0; // nnz_mat = 0;
  PetscCall(MatGetSize(A, &m_mat, &n_mat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix Row = %" PetscInt_FMT ", matrix Column = %" PetscInt_FMT "\n", m_mat, n_mat));
  MatInfo info_mat;
  PetscCall(MatGetInfo(A, MAT_GLOBAL_SUM, &info_mat));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix nz_allocated = %ld, matrix nz_used = %ld, matrix nz_unneeded = %ld\n",
                        (long)(info_mat.nz_allocated), (long)(info_mat.nz_used), (long)(info_mat.nz_unneeded)));
  MatType mtype;
  MatGetType(A, &mtype);
  PetscPrintf(PETSC_COMM_WORLD, "Matrix type = %s\n", mtype);

  const PetscInt *csr_ia = NULL;
  const PetscInt *csr_ja = NULL;
  int n_loc = 0, nnz_loc = 0;
  int loc_row_start = 0, loc_row_end = 0;
  PetscBool done = PETSC_TRUE;
  int *ja = NULL, *ia = NULL;
  double *a = NULL;

  PetscCall(MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n_loc,
                        &csr_ia, &csr_ja, &done));
  PetscCall(MatGetOwnershipRange(A, &loc_row_start, &loc_row_end));

  // n_loc = loc_row_end - loc_row_start;
  nnz_loc = csr_ia[n_loc];
  int *loc_row_idx = NULL;

#if 1
  if ((a = (double *)malloc(nnz_loc * sizeof(double))) == NULL ||
      (ja = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
      (ia = (int *)malloc((n_loc + 1) * sizeof(int))) == NULL ||
      (loc_row_idx = (int *)malloc(n_loc * sizeof(int))) == NULL)
  {
    fprintf(stderr, "Memory allocation failed! \"binary matrix\"\n");
    exit(EXIT_FAILURE);
  }
#endif

  printf("\n==== loc_row_start = %d, loc_row_end = %d\n", loc_row_start, loc_row_end);
  printf("\n==== n_loc = %d, nnz_loc = %d\n", n_loc, nnz_loc);

#if 1
  for (int index = loc_row_start; index < loc_row_end; ++index)
  {
    loc_row_idx[index - loc_row_start] = index;
  }

  for (int index = 0; index < n_loc + 1; ++index)
  {
    ia[index] = csr_ia[index];
  }

  for (int index = 0; index < nnz_loc; ++index)
  {
    ja[index] = csr_ja[index];
  }

  for (int index = 0; index < n_loc; ++index)
  {
    // matrix element
    int index_start = csr_ia[index];
    int index_end = csr_ia[index + 1];
    for (int index_j = index_start; index_j < index_end; ++index_j)
    {
      PetscScalar val_tmp;
      PetscCall(MatGetValue(A, loc_row_idx[index], csr_ja[index_j], &val_tmp));
      a[index_j] = val_tmp;
    }
  }
#endif

#if 1
  FILE *fp = fopen(output, "wb");

  if (fp == NULL)
  {
    printf("cannot open file \"%s\"\n", output);
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "%d\t%d\t%d\n", nrows, ncols, nnz_loc);
#endif

#if 1
  for (int index = 0; index < n_loc + 1; ++index)
  {
    fprintf(fp, "%d\n", ia[index]);
  }

  for (int index = 0; index < nnz_loc; ++index)
  {
    fprintf(fp, "%d\n", ja[index]);
  }

  for (int index = 0; index < nnz_loc; ++index)
  {
    fprintf(fp, "%021.16le\n", a[index]);
  }

  fclose(fp);
#endif

  PetscCall(PetscFinalize());

  return 0;
}