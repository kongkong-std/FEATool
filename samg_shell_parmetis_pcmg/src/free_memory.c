#include "../include/main.h"

int FreeSAMGMatData(MGLevel *level /*samg level data*/)
{
    if (!level)
    {
        return 0;
    }

    PetscCall(MatDestroy(&level->op_f));
    PetscCall(MatDestroy(&level->op_c));
    PetscCall(MatDestroy(&level->op_ua_p));
    PetscCall(MatDestroy(&level->op_sa_p));

    return 0;
}

int FreeQLevelK(QLevelK *q /*Q matrix data*/)
{
    if (!q)
    {
        return 0;
    }

    if (q->data_q)
    {
        free(q->data_q);
    }

    q->size_global = 0;
    q->size_local = 0;

    return 0;
}

int FreeAggData(AggData *agg /*aggregation data*/, int my_rank /*current rank*/)
{
    if (!agg)
    {
        return 0;
    }

    int np = agg->np;

    for (int index = 0; index < np; ++index)
    {
        if (agg->owner[index] == my_rank && agg->data_ghost_agg + index)
        {
            FreeGhostAggData(agg->data_ghost_agg + index);
        }
    }
    if (agg->data_ghost_agg)
    {
        free(agg->data_ghost_agg);
    }

    if (agg->fgid2part)
    {
        free(agg->fgid2part);
    }

    for (int index = 0; index < np; ++index)
    {
        if (agg->owner[index] == my_rank && agg->fine_global_nullspace[index])
        {
            free(agg->fine_global_nullspace[index]);
        }
    }
    if (agg->fine_global_nullspace)
    {
        free(agg->fine_global_nullspace);
    }

    for (int index = 0; index < np; ++index)
    {
        if (agg->owner[index] == my_rank && agg->fine_gids[index])
        {
            free(agg->fine_gids[index]);
        }
    }
    if (agg->fine_gids)
    {
        free(agg->fine_gids);
    }

    if (agg->n_fine)
    {
        free(agg->n_fine);
    }

    if (agg->nlocal_all)
    {
        free(agg->nlocal_all);
    }

    for (int index = 0; index < np; ++index)
    {
        if (agg->local_gids[index] && agg->nlocal[index] > 0)
        {
            free(agg->local_gids[index]);
        }
    }
    if (agg->local_gids)
    {
        free(agg->local_gids);
    }

    if (agg->nlocal)
    {
        free(agg->nlocal);
    }

    if (agg->owner)
    {
        free(agg->owner);
    }

    return 0;
}

int FreeGhostAggData(GhostAggData *g /*ghost data of aggregation data*/)
{
    if (!g)
    {
        return 0;
    }

    if (g->flid2fgid)
    {
        free(g->flid2fgid);
    }

    if (g->mat_t)
    {
        free(g->mat_t);
    }

    if (g->mat_q)
    {
        free(g->mat_q);
    }

    if (g->mat_r)
    {
        free(g->mat_r);
    }

    g->n_owned = 0;
    g->nghost = 0;
    g->n_total = 0;
    g->nrow = 0;
    g->ncol = 0;

    return 0;
}

int FreeNearNullSpaceLevelK(NearNullSpaceLevelK *ns /*near null space level k*/)
{
    if (!ns)
    {
        return 0;
    }

    if (ns->data_nullspace)
    {
        free(ns->data_nullspace);
    }

    ns->size_global = 0;
    ns->size_local = 0;

    return 0;
}

int FreeNearNullSpaceLevel0(NearNullSpaceLevel0 *ns /*near null space level 0*/)
{
    if (!ns)
    {
        return 0;
    }

    if (ns->data_nullspace)
    {
        free(ns->data_nullspace);
    }

    ns->size_global = 0;
    ns->size_local = 0;

    return 0;
}

int FreeMeshData(MeshData *mesh /*mesh data*/)
{
    if (!mesh)
    {
        return 0;
    }

    if (mesh->vtxdist)
    {
        free(mesh->vtxdist);
    }

    if (mesh->parts)
    {
        free(mesh->parts);
    }

    FreeMeshVtx(&mesh->data_vtx);
    FreeMeshAdj(&mesh->data_adj);

    mesh->nv = 0;
    mesh->np = 0;

    return 0;
}

int FreeMeshVtx(MeshVtx *data_vtx /*mesh vertex data*/)
{
    if (!data_vtx)
    {
        return 0;
    }

    if (data_vtx->idx)
    {
        free(data_vtx->idx);
    }

    if (data_vtx->type)
    {
        free(data_vtx->type);
    }

    if (data_vtx->data_coor)
    {
        free(data_vtx->data_coor);
    }

    data_vtx->local_nv = 0;
    data_vtx->nv = 0;

    return 0;
}

int FreeMeshAdj(MeshAdj *data_adj /*mesh adjacency data*/)
{
    if (!data_adj)
    {
        return 0;
    }

    if (data_adj->idx)
    {
        free(data_adj->idx);
    }

    if (data_adj->xadj)
    {
        free(data_adj->xadj);
    }

    if (data_adj->adjncy)
    {
        free(data_adj->adjncy);
    }

    data_adj->local_nv = 0;
    data_adj->nv = 0;

    return 0;
}
