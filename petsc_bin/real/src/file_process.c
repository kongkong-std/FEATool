#include "main.h"

void MatrixProcessSize(const char *path, Matrix *mat)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix size \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    char buffer[PETSC_MAX_PATH_LEN];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    fscanf(fp, "%d%d%d", &(mat->m), &(mat->n), &(mat->nnz));
    fclose(fp);
}

void VectorProcessSize(const char *path, Vector *vec)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - vector size \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d", &(vec->n));
    fclose(fp);
}

void MatrixProcess(const char *path, Matrix *mat, int row_start, int row_end)
{
    printf("loc_row_start = %d, loc_row_end = %d\n", row_start, row_end);
    int nnz_loc = 0;

    FILE *fp = NULL;

    // scanning file to get local nnz
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix process \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    char buffer[PETSC_MAX_PATH_LEN];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    for (int index = 0; index < mat->nnz; ++index)
    {
        int m_tmp = 0, n_tmp = 0;
        double val_tmp = 0.;
        fscanf(fp, "%d%d%lf", &m_tmp, &n_tmp, &val_tmp);
        --m_tmp;
        if (m_tmp >= row_start && m_tmp < row_end)
        {
            ++nnz_loc;
        }
    }
    printf("---- local nnz = %d\n", nnz_loc);
    fclose(fp);

    // getting local row_idx, col_idx, val
    if ((mat->row_idx = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (mat->col_idx = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (mat->val = (double *)malloc(nnz_loc * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'matrix data\'\n");
        exit(EXIT_FAILURE);
    }

    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix process \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }
    fgets(buffer, PETSC_MAX_PATH_LEN, fp);

    int loc_count = 0;
    for (int index = 0; index < mat->nnz; ++index)
    {
        int m_tmp = 0, n_tmp = 0;
        double val_tmp = 0.;
        fscanf(fp, "%d%d%lf", &m_tmp, &n_tmp, &val_tmp);
        --m_tmp;
        if (m_tmp >= row_start && m_tmp < row_end)
        {
            // base-1 to base-0
            mat->row_idx[loc_count] = m_tmp;
            mat->col_idx[loc_count] = n_tmp - 1;
            mat->val[loc_count] = val_tmp;
            ++loc_count;
        }
    }
    fclose(fp);

    // updating nnz to local nnz
    mat->nnz = nnz_loc;
}

void VectorProcess(const char *path, Vector *vec, int row_start, int row_end)
{
    printf("loc_row_start = %d, loc_row_end = %d\n", row_start, row_end);
    double *val_tmp = NULL;
    int loc_size = row_end - row_start;

    if ((val_tmp = (double *)malloc(vec->n * sizeof(double))) == NULL ||
        (vec->val = (double *)malloc(loc_size * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'vector data\'\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - vector data \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    int n_tmp = 0;
    fscanf(fp, "%d", &n_tmp);
    for (int index = 0; index < n_tmp; ++index)
    {
        fscanf(fp, "%lf", val_tmp + index);
    }

    fclose(fp);

    for (int index = row_start; index < row_end; ++index)
    {
        vec->val[index - row_start] = val_tmp[index];
    }

    // updating size to local size
    vec->n = loc_size;

    // free memory
    free(val_tmp);
}
