

#============

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine1.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
-pc_type none \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/refine1/comsol-hist-flat-none-np4.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine1.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
-pc_type hypre \
-pc_hypre_type boomeramg \
-pc_hypre_boomeramg_numfunctions 6 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/refine1/comsol-hist-flat-boomeramg-np4.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine1.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/refine1/comsol-hist-flat-rbm-np4.log



mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine1.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/refine1/comsol-hist-flat-2nd-null-np4.log

#============

