#include "../include/main.h"

int ParseConfig(MPI_Comm comm /*communicator*/,
                const char *path /*path to json file*/,
                CfgJson *config /*json config data*/)
{
    CfgFile *cfg_file = &config->cfg_file;
    CfgMG *cfg_mg = &config->cfg_mg;

    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    char *json_data = NULL;
    cJSON *root = NULL;

    if (my_rank == 0)
    {
        // root rank, file process
        FILE *fp = fopen(path, "rb");
        assert(fp);

        fseek(fp, 0, SEEK_END);
        long file_size = ftell(fp);
        fseek(fp, 0, SEEK_SET);

        json_data = (char *)malloc(file_size + 1);
        fread(json_data, 1, file_size, fp);
        json_data[file_size] = '\0';

        fclose(fp);

        // json data
        root = cJSON_Parse(json_data);
        assert(root);

        // json file
        cJSON *file_obj = cJSON_GetObjectItem(root, "file");
        if (file_obj)
        {
            cJSON *file_mat = cJSON_GetObjectItem(file_obj, "file_mat");
            cJSON *file_rhs = cJSON_GetObjectItem(file_obj, "file_rhs");
            cJSON *file_vtx = cJSON_GetObjectItem(file_obj, "file_vtx");
            cJSON *file_adj = cJSON_GetObjectItem(file_obj, "file_adj");

            if (file_mat && file_rhs && file_vtx && file_adj)
            {
                snprintf(cfg_file->file_mat, sizeof(cfg_file->file_mat), "%s", file_mat->valuestring);
                snprintf(cfg_file->file_rhs, sizeof(cfg_file->file_rhs), "%s", file_rhs->valuestring);
                snprintf(cfg_file->file_vtx, sizeof(cfg_file->file_vtx), "%s", file_vtx->valuestring);
                snprintf(cfg_file->file_adj, sizeof(cfg_file->file_adj), "%s", file_adj->valuestring);
            }
        }

        // json mla
        cJSON *mg_obj = cJSON_GetObjectItem(root, "mg");
        if (mg_obj)
        {
            cJSON *pre_smooth = cJSON_GetObjectItem(mg_obj, "pre_smooth");
            cJSON *post_smooth = cJSON_GetObjectItem(mg_obj, "post_smooth");
            cJSON *num_level = cJSON_GetObjectItem(mg_obj, "num_level");
            cJSON *num_coarse_vtx = cJSON_GetObjectItem(mg_obj, "num_coarse_vtx");
            cJSON *est_size_agg = cJSON_GetObjectItem(mg_obj, "est_size_agg");
            cJSON *ps_num_steps = cJSON_GetObjectItem(mg_obj, "ps_num_steps");
            cJSON *ps_type = cJSON_GetObjectItem(mg_obj, "ps_type");
            cJSON *ps_scale = cJSON_GetObjectItem(mg_obj, "ps_scale");

            if (pre_smooth &&
                post_smooth &&
                num_level &&
                num_coarse_vtx &&
                est_size_agg &&
                ps_num_steps &&
                ps_type &&
                ps_scale)
            {
                cfg_mg->pre_smooth = pre_smooth->valueint;
                cfg_mg->post_smooth = post_smooth->valueint;
                cfg_mg->num_level = num_level->valueint;
                cfg_mg->num_coarse_vtx = num_coarse_vtx->valueint;
                cfg_mg->est_size_agg = est_size_agg->valueint;
                cfg_mg->ps_num_steps = ps_num_steps->valueint;
                cfg_mg->ps_type = ps_type->valueint;
                cfg_mg->ps_scale = ps_type->valuedouble;
            }
        }
    }

    // broadcast data
    /*
     * config_file
     */
    (void)MPI_Bcast(config->cfg_file.file_mat, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, comm);
    (void)MPI_Bcast(config->cfg_file.file_rhs, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, comm);
    (void)MPI_Bcast(config->cfg_file.file_vtx, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, comm);
    (void)MPI_Bcast(config->cfg_file.file_adj, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, comm);

    /*
     * config_mla
     */
    (void)MPI_Bcast(&config->cfg_mg.pre_smooth, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.post_smooth, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.num_level, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.num_coarse_vtx, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.est_size_agg, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.ps_num_steps, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.ps_type, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&config->cfg_mg.ps_scale, 1, MPI_DOUBLE, 0, comm);

    // free memory
    if (my_rank == 0)
    {
        cJSON_Delete(root);
        free(json_data);
    }

    return 0;
}
