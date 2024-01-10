# SPA Applied to Harvesting Problem
The matlab files here incorporate a switch point algorithm (SPA) to numerically solve a harvesting optimal control problem. 
The optimal control problem used in these experiments came from Neubert, M.G. (2003), Marine reserves and optimal harvesting. Ecology Letters, 6: 843-849. https://doi.org/10.1046/j.1461-0248.2003.00493.x .

In order to run most of these files (except for the plotting files) one will need to download the SuiteOPT software.
This software can be accessed on https://people.clas.ufl.edu/hager/software/ . 
This software works for linux and unix systems. 
In our experiments we used the SuiteOPT version 2.0.3. You can use the 
same link but click software archive to access this older version of the software.

harvest3.m   -- Solves the harvesting problem using Total variation regularization with tuning parameter set to where a control with four switches is obtained

harvest6.m   -- Solves the harvesting problem using Total variation regularization with tuning parameter set to where a control witch six switches is obtained

SPA4.m       -- Takes the control form obtained from harvest3.m and updates the switches by using SPA

SPA4Plots.m  -- File outputs figures used in analyzing the control obtained from SPA4

SPA6.m       -- Takes the control form obtained from harvest6.m and updated the switches by using SPA

SPA6Plots.m  -- File outputs figures used in analyzing the control obtained from SPA6

SPA6bang.m   -- Takes a bang-bang control form with 6 switches and uses SPA to improve its switches

SPA6Bang


