mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/curved-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-curved-rbm-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/curved-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-curved-rbm-np4.log

mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/curved-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-curved-2nd-null-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/curved-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-curved-2nd-null-np4.log

#============

mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-flat-rbm-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-flat-rbm-np4.log

mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-flat-2nd-null-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/flat-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-flat-2nd-null-np4.log

#============

mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/inclined-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-inclined-rbm-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/inclined-shell/refine0.json \
-order_rbm 1 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-inclined-rbm-np4.log

mpirun -np 1 ./app_petsc_exe \
-config ../input/config/comsol/inclined-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-inclined-2nd-null-np1.log

mpirun -np 4 ./app_petsc_exe \
-config ../input/config/comsol/inclined-shell/refine0.json \
-order_rbm 2 \
-angle_type 1 \
-ksp_type fcg \
-ksp_max_it 2000 \
-ksp_monitor_true_residual \
-ksp_norm_type unpreconditioned \
-ksp_rtol 1e-8 \
-mesh_type 1 \
2>&1 | tee /home/kongkong/Documents/DocFEMTool/FEATool/ml_parmetis_pcmg_linsys/exe/log-comsol-hist/sa/comsol-hist-inclined-2nd-null-np4.log

#============
