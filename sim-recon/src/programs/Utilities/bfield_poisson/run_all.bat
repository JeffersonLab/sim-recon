
START C:\LANL\AUTOMESH.EXE gluex_sol_04.am

@echo Please hit return *AFTER* the AUTOMESH program finishes
@echo (after the Automesh window closes)
@PAUSE
START C:\LANL\POISSON.EXE GLUEX_SOL_04.t35

@echo Please hit return *AFTER* the POISSON program finishes (may be several minutes)
@echo (after the Poisson window closes)
@PAUSE
START C:\LANL\SF7.EXE gluex_sol_04.IN7
