mpirun -np 4 ./app_petsc_exe \
-config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type cg \
-ksp_max_it 20 \
-ksp_rtol 1e-8 \
-ksp_norm_type unpreconditioned \
-ksp_monitor_true_residual \
-def_pc_parmetis_mla \
2>&1 | tee hist-pcshell.log

#-def_pc_parmetis_mla \