mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/curved-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type cg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee comsol-hist-curved-kirchhoff.log
