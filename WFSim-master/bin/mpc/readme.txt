files from git repo:
---
- MPCinit.m: initializes controller parameters before running MPC controller
- NMPCcontroller.m: calculates optimal MPC CT_prime values 
    based on current effective wind speeds and thrust coeffecients 
- SMPCcontroller.m: calculates optimal MPC CT_prime values based on current
    effective wind speeds and thrust coeffecients
- matsys.m: sets up the stacked matrices AA, BBt and CC 
- wfmodel.m: wfmodel builds input and system matrices for wind farm model

qLMPC reformulated as QP
--
- MPCinit0.m: initializes controller parameters before running closed-loop
    simulation
- matsysLSHg.m: sets up the stacked matrices L and S as well as Hessian H and 
    gradient g    
- myquadprog.m:  myquadprog Quadratic programming based on QUADPROG with 
    flexibility removed for speed-up 
- qLMPCcontroller.m: calculates optimal MPC CT_prime values based on current
    effective wind speeds and thrust coeffecients via QP  

velocity qLMPC:
---
ABqL_WF.m -       
ABqLd_WF.m -     
HSqLd_WF.m -


test functions
---
qLMPCtrlTest.m: compares NMPCcontroller.m vs qLMPCcontroller.m
qLMPCtestHg.m:    
  
      