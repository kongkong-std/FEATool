mpirun -np 4 ./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee  log-hist/kirchhoff/refine0/hist.log

mpirun -np 4 ./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee  log-hist/kirchhoff/refine3/hist.log

mpirun -np 4 ./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee  log-hist/kirchhoff/refine2/hist.log

mpirun -np 4 ./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 2>&1 | tee  log-hist/kirchhoff/refine1/hist.log
