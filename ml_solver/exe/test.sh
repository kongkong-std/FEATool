./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 1 | tee log/refine0/rbm1-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 1 | tee log/refine1/rbm1-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 1 | tee log/refine2/rbm1-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 1 | tee log/refine3/rbm1-multilevel.log

./app_petsc_exe -config ../input/config/mla_v_cycle-refine0.json -order_rbm 2 | tee log/refine0/rbm2-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine1.json -order_rbm 2 | tee log/refine1/rbm2-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine2.json -order_rbm 2 | tee log/refine2/rbm2-multilevel.log
./app_petsc_exe -config ../input/config/mla_v_cycle-refine3.json -order_rbm 2 | tee log/refine3/rbm2-multilevel.log
