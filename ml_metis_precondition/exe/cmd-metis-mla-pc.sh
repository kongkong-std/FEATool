./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm1-refine0.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm2-refine0.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm1-refine3.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm2-refine3.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm1-refine2.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm2-refine2.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm1-refine1.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-def_pc_metis_mla \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee hist-metis-mlasolver-rbm2-refine1.log
