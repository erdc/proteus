-rans2p_ksp_type gmres -rans2p_pc_type asm -rans2p_pc_asm_type basic -rans2p_ksp_max_it 2000
-log_summary
-clsvof_ksp_type fgmres -clsvof_pc_type hypre -clsvof_pc_hypre_type boomeramg -clsvof_ksp_max_it 2000
-phi_ksp_type fgmres -phi_pc_type  hypre -phi_pc_hypre_type boomeramg -phi_ksp_max_it 2000
-pressure_ksp_type fgmres -pressure_pc_type   hypre -pressure_pc_hypre_type    boomeramg -pressure_ksp_gmres_restart 300 -pressure_ksp_knoll -pressure_ksp_max_it 2000
-pinit_ksp_type fgmres -pinit_pc_type   hypre -pinit_pc_hypre_type    boomeramg -pinit_ksp_gmres_restart 300 -pinit_ksp_knoll -pinit_ksp_max_it 2000
####
#-rans3p_ksp_type fgmres -rans3p_pc_type hypre -rans3p_pc_hypre_type boomeramg -rans3p_ksp_max_it 2000
#-rans3p_ksp_type fgmres -rans3p_pc_type hypre -rans3p_pc_hypre_type pilut -rans3p_ksp_max_it 2000
#-rans3p_ksp_type fgmres -rans3p_pc_type asm -rans3p_pc_asm_type basic -rans3p_ksp_max_it 2000
#-rans3p_ksp_type fgmres -rans3p_pc_type jacobi -rans3p_ksp_max_it 2000
#-rans3p_ksp_type bicg -rans3p_ksp_max_it 2000
-rans3p_ksp_type fgmres -rans3p_pc_type asm -rans3p_pc_asm_type basic -rans3p_ksp_max_it 2000
-rans3p_ksp_gmres_modifiedgramschmidt -rans3p_ksp_gmres_restart 300 -rans3p_sub_ksp_type preonly -rans3p_sub_pc_factor_mat_solver_package superlu -rans3p_ksp_knoll -rans3p_sub_pc_type lu 
