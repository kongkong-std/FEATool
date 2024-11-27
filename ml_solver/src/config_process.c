#include "../include/main.h"

// 读取和解析 JSON 配置文件
int ConfigParse(const char *filename, ConfigFile *file_config, ConfigMLA *mla_config)
{
    // 读取 JSON 文件内容
    FILE *f = fopen(filename, "rb");
    if (!f)
    {
        fprintf(stderr, "Unable to open file %s\n", filename);
        return -1;
    }

    fseek(f, 0, SEEK_END);
    long file_size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *json_data = (char *)malloc(file_size + 1);
    fread(json_data, 1, file_size, f);
    json_data[file_size] = '\0';
    fclose(f);

    // 解析 JSON 数据
    cJSON *root = cJSON_Parse(json_data);
    if (root == NULL)
    {
        fprintf(stderr, "Error parsing JSON\n");
        free(json_data);
        return -1;
    }

    // 解析 file 配置
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

    // 解析 mla 配置
    cJSON *mla_obj = cJSON_GetObjectItem(root, "mla");
    if (mla_obj)
    {
        cJSON *pre_smooth_v = cJSON_GetObjectItem(mla_obj, "pre_smooth_v");
        cJSON *post_smooth_v = cJSON_GetObjectItem(mla_obj, "post_smooth_v");
        cJSON *mla_max_it = cJSON_GetObjectItem(mla_obj, "mla_max_it");
        cJSON *mla_rtol = cJSON_GetObjectItem(mla_obj, "mla_rtol");
        cJSON *mla_level = cJSON_GetObjectItem(mla_obj, "mla_level");
        cJSON *mla_phase = cJSON_GetObjectItem(mla_obj, "mla_phase");

        if (pre_smooth_v && post_smooth_v && mla_max_it && mla_rtol && mla_level && mla_phase)
        {
            mla_config->pre_smooth_v = pre_smooth_v->valueint;
            mla_config->post_smooth_v = post_smooth_v->valueint;
            mla_config->mla_max_it = mla_max_it->valueint;
            mla_config->mla_rtol = mla_rtol->valuedouble;
            mla_config->mla_level = mla_level->valueint;
            mla_config->mla_phase = mla_phase->valueint;
        }
    }

    // 释放分配的内存
    cJSON_Delete(root);
    free(json_data);

    return 0;
}