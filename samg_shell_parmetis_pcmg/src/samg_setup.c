#include "../include/main.h"

int SAMGApplyProlongationSmoother(PSmoother *p_s /*prolongation operator smoother*/,
                                  Mat *p_sa /*smoothed prolongation operator*/,
                                  Mat *p_ua /*unsmoothed prolongation operator*/)
{
    return 0;
}

int SAMGSetupPhase(SAMGCtx *samg_ctx /*samg context data*/)
{
    int cfg_mg_num_level = samg_ctx->data_cfg.cfg_mg.num_level;
    int cfg_mg_pre_smooth = samg_ctx->data_cfg.cfg_mg.pre_smooth;
    int cfg_mg_post_smooth = samg_ctx->data_cfg.cfg_mg.post_smooth;
    int cfg_mg_num_coarse_vtx = samg_ctx->data_cfg.cfg_mg.num_coarse_vtx;
    int cfg_mg_est_size_agg = samg_ctx->data_cfg.cfg_mg.est_size_agg;
    int cfg_mg_ps_num_steps = samg_ctx->data_cfg.cfg_mg.ps_num_steps;
    int cfg_mg_ps_type = samg_ctx->data_cfg.cfg_mg.ps_type;
    double cfg_mg_ps_scale = samg_ctx->data_cfg.cfg_mg.ps_scale;

    samg_ctx->levels = (MGLevel *)malloc(cfg_mg_num_level * sizeof(MGLevel));
    assert(samg_ctx->levels);

    int cnt_level = 0;
    PetscCall(MatDuplicate(samg_ctx->mysolver.solver_a, MAT_COPY_VALUES, &samg_ctx->levels[cnt_level].op_f));
    while (cnt_level < cfg_mg_num_level)
    {
        ++cnt_level;
    }

    samg_ctx->num_level = cnt_level;

    return 0;
}
