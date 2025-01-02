module purge
module load compilers/gcc/10.3.0 generators/make/4.3/gcc-10.3.0 mpi/mpich/3.2.1/gcc-10.3.0 tools/cJSON/1.7.18

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type none | tee log/fgmres/refine0/no-pre/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type jacobi | tee log/fgmres/refine0/jacobi/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type sor | tee log/fgmres/refine0/sor/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type ilu | tee log/fgmres/refine0/ilu/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 1 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine0/boomeramg-num1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 6 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine0/boomeramg-num6/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine0/mla-rbm1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 2 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine0/mla-rbm2/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type none | tee log/fgmres/refine3/no-pre/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type jacobi | tee log/fgmres/refine3/jacobi/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type sor | tee log/fgmres/refine3/sor/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type ilu | tee log/fgmres/refine3/ilu/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 1 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine3/boomeramg-num1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 6 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine3/boomeramg-num6/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine3/mla-rbm1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 2 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine3/mla-rbm2/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type none | tee log/fgmres/refine2/no-pre/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type jacobi | tee log/fgmres/refine2/jacobi/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type sor | tee log/fgmres/refine2/sor/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type ilu | tee log/fgmres/refine2/ilu/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 1 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine2/boomeramg-num1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 6 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine2/boomeramg-num6/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine2/mla-rbm1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 2 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine2/mla-rbm2/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type none | tee log/fgmres/refine1/no-pre/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type jacobi | tee log/fgmres/refine1/jacobi/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type sor | tee log/fgmres/refine1/sor/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type ilu | tee log/fgmres/refine1/ilu/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 1 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine1/boomeramg-num1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-pc_type hypre -pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 6 -pc_hypre_boomeramg_max_levels 2 | tee log/fgmres/refine1/boomeramg-num6/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine1/mla-rbm1/result.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 2 \
-ksp_type fgmres \
-ksp_gmres_restart 50 \
-ksp_max_it 2000 \
-ksp_rtol 1e-8 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-def_pc_mla | tee log/fgmres/refine1/mla-rbm2/result.log
