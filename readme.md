Directories
- WFSimEnvironment: Contains modified WFsim environment for farm simulation
  open and closed loop
- KoopmanIODMD: Contains modified Koopman sysId

First step: To generate OpenLoop Data for testing (in folder WFSimEnvironment/Code/)
- run WFSim_demo with controller settings from controlDefinitions and controller = 0
-- saves data in format expected by MainWFSim into KoopmanIODMD/DataT2OLWFSim 
 Ct1, Ct2,FT1,FT2,PT1,PT2,PtotalGreedy,Ur1,Ur2,phi1,phi2 ,u , v,p  

Second step: To generate Koopman Sys Id for windfarm controller 
- run KoopmanIODMD/Code/MainWfSim:
-- saves Koopman SysID model in subfolders of  in KoopmanIODMD/DataInOutWfSim 

Third step: Evaluate quality of Koopman Sys ID in WFsim in closed loop with MPC
- run code/WFSim_demo with wind speed used for data generation and controller = 2

Add YALMIP to run code- https://github.com/yalmip/YALMIP 

References: 
1-WFSim_demo is the WIND FARM SIMULATOR (WFSIM) by S. Boersma and B. Doekemeijer
 from https://github.com/TUDelft-DataDrivenControl/WFSimhttps://github.com/TUDelft-DataDrivenControl/WFSim
2- Cassamo, N., & van Wingerden, J. W. (2021, May). Model Predictive Control for Wake Redirection in Wind Farms: a Koopman Dynamic Mode Decomposition Approach. 
In 2021 American Control Conference (ACC) (pp. 1776-1782). IEEE.
https://github.com/nassircassamo/IODMD_SOWFA

try
