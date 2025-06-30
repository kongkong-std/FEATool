#include <libMLA.h>

int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscBool path_flag;

    char path_config[PETSC_MAX_PATH_LEN];
    PetscInt order_rbm = 0, angle_type = 0;

    PetscCall(PetscOptionsGetString(NULL, NULL, "-config", path_config, sizeof(path_config), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(comm, "Config file: %s\n", path_config));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-order_rbm", &order_rbm, &path_flag));
    if (path_flag)
    {
        // parameter order_rbm
        /*
         * order_rbm = 1, classical rbm type prolongation operator
         * order_rbm = 2, krichhoff type prolongation operator
         */
        PetscCall(PetscPrintf(comm, "order_rbm = %" PetscInt_FMT "\n", order_rbm));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-angle_type", &angle_type, &path_flag));
    if (path_flag)
    {
        // parameter angle_type
        /*
         * angle_type = 0, rotational angle with axis
         * angle_type = 1, normal rotational angle with axis
         */
        PetscCall(PetscPrintf(comm, "angle_type = %" PetscInt_FMT "\n", angle_type));
    }

    // user-defined precondition
    PetscBool def_pc_parmetis_mla = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_parmetis_mla", &def_pc_parmetis_mla, NULL));

    // mla_ctx
    MLAContext mla_ctx;
    mla_ctx.setup = 0;
    mla_ctx.path_config = path_config;
    mla_ctx.order_rbm = order_rbm;
    mla_ctx.angle_type = angle_type;

    int status_config_json = ConfigParse(comm, mla_ctx.path_config, &(mla_ctx.config));
    if (status_config_json != 0)
    {
        exit(EXIT_FAILURE);
    }
    mla_ctx.metis_mla = (MetisMLAGraph *)malloc(mla_ctx.config.mla_config.mla_level *
                                                sizeof(MetisMLAGraph));
    assert(mla_ctx.metis_mla);
#if 0
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            puts("\n==== mla_ctx information ====\n");
            printf("in rank %d/%d\n", my_rank, nprocs);
            printf("setup = %d\n", mla_ctx.setup);
            printf("path_config = %s\n", mla_ctx.path_config);
            printf("order_rbm = %d\n", mla_ctx.order_rbm);
            printf("angle_type = %d\n", mla_ctx.angle_type);

            puts("\n--------\n");

            printf("file_mat = %s\n", mla_ctx.config.file_config.file_mat);
            printf("file_rhs = %s\n", mla_ctx.config.file_config.file_rhs);
            printf("file_mesh = %s\n", mla_ctx.config.file_config.file_mesh);
            printf("pre_smooth_v = %d\n", mla_ctx.config.mla_config.pre_smooth_v);
            printf("post_smooth_v = %d\n", mla_ctx.config.mla_config.post_smooth_v);
            printf("mla_max_it = %d\n", mla_ctx.config.mla_config.mla_max_it);
            printf("mla_rtol = %021.1le\n", mla_ctx.config.mla_config.mla_rtol);
            printf("mla_level = %d\n", mla_ctx.config.mla_config.mla_level);
            printf("mla_phase = %d\n", mla_ctx.config.mla_config.mla_phase);
            printf("coarse_restart = %d\n", mla_ctx.config.mla_config.coarse_restart);
            printf("label_bound = %d\n", mla_ctx.config.mesh_label_config.label_bound);
            printf("label_omega = %d\n", mla_ctx.config.mesh_label_config.label_omega);

            puts("\n--------\n");
        }
    }
#endif // mla_ctx data information

    /*
     * step 1, mesh file process
     */
    DataMesh mesh_data;

    int *node_vtxdist = (int *)calloc(nprocs + 1, sizeof(int));
    assert(node_vtxdist);

    if (my_rank == 0)
    {
        FileProcessMesh(mla_ctx.config.file_config.file_mesh, &mesh_data);
#if 0
        puts("\n>>>> information of mesh >>>>\n");
        printf("nn of mesh: %d\n", mesh_data.nn);
        for (int index = 0; index < mesh_data.nn; ++index)
        {
            printf("node %d:\t", index);
            for (int index_i = 0; index_i < mesh_data.dim; ++index_i)
            {
                printf("%021.16le\t", mesh_data.coordinates[mesh_data.dim * index + index_i]);
            }
            putchar('\n');
        }

        printf("\nne of shell: %d\n", mesh_data.ne_shell);
        for (int index = 0; index < mesh_data.ne_shell; ++index)
        {
            printf("element %d:\t", index);
            for (int index_i = 0; index_i < mesh_data.ele_shell[index].num_ele_node; ++index_i)
            {
                printf("%d\t", mesh_data.ele_shell[index].ele_node[index_i]);
            }
            putchar('\n');
        }

        puts("\n>>>>>>>>\n\n");
#endif // mesh data

        // csr format graph
        GlobalGraphCSRAdjGenerator(&mesh_data, &mla_ctx.data_mesh);
#if 0
        puts("\n>>>> information of graph >>>>\n");
        for (int index = 0; index < mla_ctx.data_mesh->nn; ++index)
        {
            printf("node %d:\t", index);
            for (int index_i = 0; index_i < mla_ctx.data_mesh->dim; ++index_i)
            {
                printf("%021.16le\t",
                       mla_ctx.data_mesh->coordinates[mla_ctx.data_mesh->dim * index + index_i]);
            }
            putchar('\n');
        }

        puts("\n\ngraph_data xadj value:");
        for (int index = 0; index < mla_ctx.data_mesh->nn + 1; ++index)
        {
            printf("%" PRIDX "\t", mla_ctx.data_mesh->xadj[index]);
        }
        putchar('\n');

        puts("\n\ngraph_data adjncy value:");
        for (int index = 0; index < mla_ctx.data_mesh->xadj[mla_ctx.data_mesh->nn]; ++index)
        {
            printf("%" PRIDX "\t", mla_ctx.data_mesh->adjncy[index]);
        }
        putchar('\n');

        puts("\n>>>>>>>>\n\n");
#endif // graph csr data

        int base_num = mla_ctx.data_mesh->nn / nprocs;
        int remainder_num = mla_ctx.data_mesh->nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            node_vtxdist[index_p + 1] = node_vtxdist[index_p] + count;
        }
    }

    (void)MPI_Bcast(node_vtxdist,
                    (nprocs + 1) * sizeof(int),
                    MPI_BYTE, 0, comm);
#if 0
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        (void)MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("in rank %d/%d\n", my_rank, nprocs);
            for (int index = 0; index < nprocs + 1; ++index)
            {
                printf("%d\t", node_vtxdist[index]);
            }
            puts("\n\n========\n\n");
        }
    }
#endif // print vtxdist data

    /*
     * step 2, linear system assembling
     */
    MySolver mysolver;
    SolverPetscInitialize(mla_ctx.config.file_config.file_mat,
                          mla_ctx.config.file_config.file_rhs,
                          node_vtxdist,
                          &mysolver);

    PetscInt local_vec_size = 0;
    PetscCall(VecGetLocalSize(mysolver.solver_b, &local_vec_size));
    printf(">>>>local_vec_size = %d\n", local_vec_size);

    PetscInt local_mat_m = 0, local_mat_n = 0;
    PetscCall(MatGetLocalSize(mysolver.solver_a, &local_mat_m, &local_mat_n));
    printf(">>>>~~~~ local_mat_m = %d, local_mat_n = %d\n", local_mat_m, local_mat_n);

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx.mysolver.ksp)));
    PetscCall(DeepCopyMLAContextMySolver(&(mla_ctx.mysolver), &mysolver));

    SolverPetscResidualCheck(&mla_ctx.mysolver);

#if 0
    // file IO
    FILE *fp = fopen("mat.txt", "wb");
    assert(fp);

    /*
     * matrix
     */
    PetscInt nrows_mat, ncols_mat, nnzs_mat;
    PetscCall(MatGetSize(mysolver.solver_a, &nrows_mat, &ncols_mat));

    fprintf(fp, "%d\t%d\n", nrows_mat, ncols_mat);

    PetscInt *row_ptr = (PetscInt *)calloc(nrows_mat + 1, sizeof(PetscInt));
    assert(row_ptr);

    for (int index = 0; index < nrows_mat; ++index)
    {
        PetscInt nnz_row = 0;
        PetscCall(MatGetRow(mysolver.solver_a, index, &nnz_row,
                            NULL, NULL));
        row_ptr[index + 1] = row_ptr[index] + nnz_row;
    }

    nnzs_mat = row_ptr[nrows_mat];

    fprintf(fp, "%d\n", nnzs_mat);

    for (int index = 0; index < nrows_mat + 1; ++index)
    {
        fprintf(fp, "%d\n", row_ptr[index]);
    }

    PetscInt *col_idx = (PetscInt *)malloc(nnzs_mat * sizeof(PetscInt));
    PetscScalar *val_mat = (PetscScalar *)malloc(nnzs_mat * sizeof(PetscScalar));
    assert(col_idx && val_mat);

    int tmp_cnt = 0;
    for (int index = 0; index < nrows_mat; ++index)
    {
        PetscInt nnz_row = 0;
        PetscInt *cols_row;
        PetscScalar *vals_row;
        PetscCall(MatGetRow(mysolver.solver_a, index, &nnz_row,
                            (const PetscInt **)&cols_row,
                            (const PetscScalar **)&vals_row));

        for (int index_i = 0; index_i < nnz_row; ++index_i)
        {
            col_idx[tmp_cnt] = cols_row[index_i];
            val_mat[tmp_cnt] = vals_row[index_i];
            ++tmp_cnt;
        }

        PetscCall(MatRestoreRow(mysolver.solver_a, index, &nnz_row,
                                (const PetscInt **)&cols_row,
                                (const PetscScalar **)&vals_row));
    }

    printf("\n\n>>>>tmp_cnt = %d\n\n", tmp_cnt);

    for (int index = 0; index < nnzs_mat; ++index)
    {
        fprintf(fp, "%d\n", col_idx[index]);
    }

    for (int index = 0; index < nnzs_mat; ++index)
    {
        fprintf(fp, "%021.16le\n", val_mat[index]);
    }

    fclose(fp);

    fp = fopen("rhs.txt", "wb");
    assert(fp);

    /*
     * rhs
     */
    PetscInt nrows_rhs;
    PetscCall(VecGetSize(mysolver.solver_b, &nrows_rhs));

    fprintf(fp, "%d\n", nrows_rhs);

    const PetscScalar *rhs_array;
    PetscCall(VecGetArrayRead(mysolver.solver_b, &rhs_array));
    for (PetscInt i = 0; i < nrows_rhs; ++i)
    {
        fprintf(fp, "%021.16le\n", rhs_array[i]);
    }
    PetscCall(VecRestoreArrayRead(mysolver.solver_b, &rhs_array));

    fclose(fp);
#endif // linsys file IO, serial version

    /*
     * test
     * calling ParMetisMLASolver
     */
    PetscCall(ParMetisMLASolver(&mla_ctx, 2));

    /*
     * linear system solver
     */
    PetscCall(KSPSetOperators(mysolver.ksp, mysolver.solver_a, mysolver.solver_a));
    PetscCall(KSPGetPC(mysolver.ksp, &(mysolver.pc)));
    if (def_pc_parmetis_mla)
    {
        PetscCall(PCSetType(mysolver.pc, PCSHELL));
        PetscCall(PCShellSetContext(mysolver.pc, &mla_ctx));
        PetscCall(PCShellSetSetUp(mysolver.pc, ParMetisMLAShellPCSetup));
        PetscCall(PCShellSetApply(mysolver.pc, ParMetisMLAShellPCApply));
    }

    PetscCall(KSPSetFromOptions(mysolver.ksp));

    PetscLogDouble time1, time2;
    PetscCall(PetscTime(&time1));
    PetscCall(KSPSolve(mysolver.ksp, mysolver.solver_b, mysolver.solver_x));
    PetscCall(PetscTime(&time2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> ksp solve time: %g (s)\n", time2 - time1));

    SolverPetscResidualCheck(&mysolver);

    // free memory
    free(mla_ctx.metis_mla);
    if (my_rank == 0)
    {
        // graph data
        free(mla_ctx.data_mesh->adjncy);
        free(mla_ctx.data_mesh->xadj);
        free(mla_ctx.data_mesh->coordinates);
        free(mla_ctx.data_mesh);

        // mesh data
        free(mesh_data.coordinates);
        for (int index = 0; index < mesh_data.ne_shell; ++index)
        {
            free(mesh_data.ele_shell[index].ele_node);
        }
        free(mesh_data.ele_shell);
    }

    PetscCall(MatDestroy(&(mla_ctx.mysolver.solver_a)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_b)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_x)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_r)));
    PetscCall(KSPDestroy(&mla_ctx.mysolver.ksp));
    PetscCall(KSPDestroy(&(mysolver.ksp)));
    PetscCall(MatDestroy(&(mysolver.solver_a)));
    PetscCall(VecDestroy(&(mysolver.solver_b)));
    PetscCall(VecDestroy(&(mysolver.solver_x)));
    PetscCall(VecDestroy(&(mysolver.solver_r)));

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}



#if 0
int DeepCopyMLAContextMySolver(MySolver *dst, MySolver *src);

int main(int argc, char **argv)
{

#if 0
    puts("\n==== test metis function ====");
    TestMetisFunctionGmsh(*(mla_ctx.data_gmsh));
#endif // test metis function

    /*
     * step 3, calling MetisMLASolver()
     */
    // PetscCall(MetisMLASolver(&mla_ctx, 2));

#if 1
    // pc shell
    PetscCall(KSPSetOperators(mysolver.ksp, mysolver.solver_a, mysolver.solver_a));
    PetscCall(KSPGetPC(mysolver.ksp, &(mysolver.pc)));
    if (def_pc_mla)
    {
        PetscCall(PCSetType(mysolver.pc, PCSHELL));
        PetscCall(PCShellSetContext(mysolver.pc, &mla_ctx));
        PetscCall(PCShellSetSetUp(mysolver.pc, MLAShellPCSetup));
        PetscCall(PCShellSetApply(mysolver.pc, MLAShellPCApply));
    }

    if (def_pc_metis_mla)
    {
        PetscCall(PCSetType(mysolver.pc, PCSHELL));
        PetscCall(PCShellSetContext(mysolver.pc, &mla_ctx));
        PetscCall(PCShellSetSetUp(mysolver.pc, MetisMLAShellPCSetup));
        PetscCall(PCShellSetApply(mysolver.pc, MetisMLAShellPCApply));
    }

    PetscCall(KSPSetFromOptions(mysolver.ksp));

    PetscLogDouble time1, time2;
    PetscCall(PetscTime(&time1));
    PetscCall(KSPSolve(mysolver.ksp, mysolver.solver_b, mysolver.solver_x));
    PetscCall(PetscTime(&time2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> ksp solve time: %g (s)\n", time2 - time1));
#endif // pcshell

   

    // check, computing residual
    SolverPetscResidualCheck(&mysolver);

    // free memory
    if (def_pc_mla)
    {
        for (int index = 0; index < mla_ctx.num_level; ++index)
        {
            free(mla_ctx.mla[index].coarse_node);
            PetscCall(MatDestroy(&(mla_ctx.mla[index].prolongation)));
            PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_coarse)));
            PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_fine)));
            PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_coarse)));
            PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_presmooth)));
            PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_postsmooth)));
            PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_coarse)));
            PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_presmooth)));
            PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_postsmooth)));
            ClearMeshGraph(mla_ctx.mla[index].aggregation);
            ClearMeshGraph(mla_ctx.mla[index].fine);
            ClearMeshGraph(mla_ctx.mla[index].coarse);
            ClearMeshGraph(mla_ctx.mla[index].mesh_tmp);
        }
        free(mla_ctx.mla);
    }

    if (def_pc_metis_mla)
    {
        for (int index = 0; index < mla_ctx.num_level; ++index)
        {
            free(mla_ctx.metis_mla[index].fine->coordinates);
            free(mla_ctx.metis_mla[index].fine->eptr_in);
            free(mla_ctx.metis_mla[index].fine->eind_in);
            free(mla_ctx.metis_mla[index].fine->npart_in);
            free(mla_ctx.metis_mla[index].fine);

            free(mla_ctx.metis_mla[index].coarse->coordinates);
            free(mla_ctx.metis_mla[index].coarse->eptr_in);
            free(mla_ctx.metis_mla[index].coarse->eind_in);
            free(mla_ctx.metis_mla[index].coarse->npart_in);
            free(mla_ctx.metis_mla[index].coarse);

            PetscCall(MatDestroy(&(mla_ctx.metis_mla[index].operator_coarse)));
            PetscCall(MatDestroy(&(mla_ctx.metis_mla[index].operator_fine)));
            PetscCall(MatDestroy(&(mla_ctx.metis_mla[index].prolongation)));
        }
    }

    ClearList(&data_list_ele_omega);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_node);
    ClearList(&data_list_phy_tag);
    ClearMeshGraph(mla_ctx.graph);

    free(mla_ctx.data_gmsh->eptr_bd);
    free(mla_ctx.data_gmsh->eptr_in);
    free(mla_ctx.data_gmsh->eind_bd);
    free(mla_ctx.data_gmsh->eind_in);
    free(mla_ctx.data_gmsh->epart_in);
    free(mla_ctx.data_gmsh->npart_in);
    free(mla_ctx.data_gmsh->coordinates);
    free(mla_ctx.data_gmsh);
    free(mla_ctx.metis_mla);

    PetscCall(MatDestroy(&(mla_ctx.mysolver.solver_a)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_b)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_x)));
    PetscCall(VecDestroy(&(mla_ctx.mysolver.solver_r)));
    PetscCall(KSPDestroy(&(mysolver.ksp)));
    PetscCall(MatDestroy(&(mysolver.solver_a)));
    PetscCall(VecDestroy(&(mysolver.solver_b)));
    PetscCall(VecDestroy(&(mysolver.solver_x)));
    PetscCall(VecDestroy(&(mysolver.solver_r)));
#if 0
    free(mla_ctx.data_gmsh->coordinates);
    free(mla_ctx.data_gmsh->eptr_bd);
    free(mla_ctx.data_gmsh->eind_bd);
    free(mla_ctx.data_gmsh->eptr_in);
    free(mla_ctx.data_gmsh->eind_in);
    free(mla_ctx.data_gmsh->epart_in);
    free(mla_ctx.data_gmsh->npart_in);
    free(mla_ctx.data_gmsh);
    for (int index = 0; index < mla_ctx.num_level; ++index)
    {
        free(mla_ctx.mla[index].coarse_node);
        PetscCall(MatDestroy(&(mla_ctx.mla[index].prolongation)));
        PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_coarse)));
        PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_fine)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_coarse)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_presmooth)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_postsmooth)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_coarse)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_presmooth)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_postsmooth)));
        ClearMeshGraph(mla_ctx.mla[index].aggregation);
        ClearMeshGraph(mla_ctx.mla[index].fine);
        ClearMeshGraph(mla_ctx.mla[index].coarse);
        ClearMeshGraph(mla_ctx.mla[index].mesh_tmp);
    }
    free(mla_ctx.mla);
    ClearList(&data_list_ele_omega);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_node);
    ClearList(&data_list_phy_tag);
    ClearMeshGraph(mla_ctx.graph);
    PetscCall(MatDestroy(&mysolver.solver_a));
    PetscCall(VecDestroy(&mysolver.solver_b));
    PetscCall(VecDestroy(&mysolver.solver_x));
    PetscCall(VecDestroy(&mysolver.solver_r));
    PetscCall(KSPDestroy(&mysolver.ksp));
    // PetscCall(PCDestroy(&mysolver.pc));
#endif // free memory

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}
#endif // serial version

// command line
/* >>>> for mla solver case
 * ./app_petsc_exe -config </path/to/config/file>
 *     -order_rbm <int, 1st order, 2nd order>
 */

/*
 * >>>> for mla preconditioning case
 * ./app_petsc_exe -config </path/to/config/file>
 *     -order_rbm <int, 1st order, 2nd order>
 *     -ksp_type fgmres
 *     -def_pc_mla
 *     -ksp_max_it 2000
 *     -ksp_rtol 1e-8
 *     -ksp_gmres_restart 50
 *     -ksp_monitor_true_residual
 *     -ksp_norm_type unpreconditioned
 *     -pc_type hypre
 *     -pc_hypre_type boomeramg
 *     -pc_hypre_boomeramg_numfunctions
 */
