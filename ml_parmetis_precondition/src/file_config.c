#include "../include/main.h"

int ConfigParse(MPI_Comm comm /*communicator*/,
                const char *path /*path to json file*/,
                ConfigJSON *config /*json config data*/)
{
    ConfigFile *file_config = &config->file_config;
    ConfigMLA *mla_config = &config->mla_config;
    ConfigMeshLabel *mesh_label_config = &config->mesh_label_config;

    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    if (my_rank == 0)
    {
        // root rank, file process
        FILE *fp = fopen(path, "rb");
        assert(fp);

        fseek(fp, 0, SEEK_END);
        long file_size = ftell(fp);
        fseek(fp, 0, SEEK_SET);

        char *json_data = (char *)malloc(file_size + 1);
        fread(json_data, 1, file_size, fp);
        json_data[file_size] = '\0';

        fclose(fp);

        // json data
        cJSON *root = cJSON_Parse(json_data);
        assert(root);

        // json file
        cJSON *file_obj = cJSON_GetObjectItem(root, "file");
        if (file_obj)
        {
            cJSON *file_mat = cJSON_GetObjectItem(file_obj, "file_mat");
            cJSON *file_rhs = cJSON_GetObjectItem(file_obj, "file_rhs");
            cJSON *file_mesh = cJSON_GetObjectItem(file_obj, "file_mesh");

            if (file_mat && file_rhs && file_mesh)
            {
                snprintf(file_config->file_mat, sizeof(file_config->file_mat), "%s", file_mat->valuestring);
                snprintf(file_config->file_rhs, sizeof(file_config->file_rhs), "%s", file_rhs->valuestring);
                snprintf(file_config->file_mesh, sizeof(file_config->file_mesh), "%s", file_mesh->valuestring);
            }
        }

        // json mla
        cJSON *mla_obj = cJSON_GetObjectItem(root, "mla");
        if (mla_obj)
        {
            cJSON *pre_smooth_v = cJSON_GetObjectItem(mla_obj, "pre_smooth_v");
            cJSON *post_smooth_v = cJSON_GetObjectItem(mla_obj, "post_smooth_v");
            cJSON *mla_max_it = cJSON_GetObjectItem(mla_obj, "mla_max_it");
            cJSON *mla_rtol = cJSON_GetObjectItem(mla_obj, "mla_rtol");
            cJSON *mla_level = cJSON_GetObjectItem(mla_obj, "mla_level");
            cJSON *mla_phase = cJSON_GetObjectItem(mla_obj, "mla_phase");
            cJSON *coarse_restart = cJSON_GetObjectItem(mla_obj, "coarse_restart");

            if (pre_smooth_v &&
                post_smooth_v &&
                mla_max_it &&
                mla_rtol &&
                mla_level &&
                mla_phase &&
                coarse_restart)
            {
                mla_config->pre_smooth_v = pre_smooth_v->valueint;
                mla_config->post_smooth_v = post_smooth_v->valueint;
                mla_config->mla_max_it = mla_max_it->valueint;
                mla_config->mla_rtol = mla_rtol->valuedouble;
                mla_config->mla_level = mla_level->valueint;
                mla_config->mla_phase = mla_phase->valueint;
                mla_config->coarse_restart = coarse_restart->valueint;
            }
        }

        // json mesh label
        cJSON *mesh_label_obj = cJSON_GetObjectItem(root, "mesh_label");
        if (mesh_label_obj)
        {
            cJSON *label_bound = cJSON_GetObjectItem(mesh_label_obj, "label_bound");
            cJSON *label_omega = cJSON_GetObjectItem(mesh_label_obj, "label_omega");

            if (label_bound && label_omega)
            {
                mesh_label_config->label_bound = label_bound->valueint;
                mesh_label_config->label_omega = label_omega->valueint;
            }
        }
    }

    // free memory
    if (my_rank == 0)
    {
        cJSON_Delete(root);
        free(json_data);
    }

    return 0;
}
