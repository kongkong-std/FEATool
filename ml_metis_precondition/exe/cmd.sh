#### order_rbm = 1
#### angle_type = 1
./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8



./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8




./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


#############################
#### order_rbm = 2
#### angle_type = 1


./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8



./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8




./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8





########################################
# order_rbm = 1
# angle_type = 0

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 1 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8



./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 1 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8




./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 1 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 1 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


#############################
# order_rbm = 2
# angle_type = 0


./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json \
-order_rbm 2 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8



./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json \
-order_rbm 2 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8




./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json \
-order_rbm 2 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8


./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json \
-order_rbm 2 \
-angle_type 0 \
-ksp_type fcg \
-ksp_max_it 20 \
-ksp_rtol 1e-8
