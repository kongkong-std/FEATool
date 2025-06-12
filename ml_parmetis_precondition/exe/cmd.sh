mpirun -np 4 ./app_petsc_exe \
-config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 1 \
-angle_type 1 \
2>&1 | tee hist.log