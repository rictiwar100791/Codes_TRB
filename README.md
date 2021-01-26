# Codes_TRB
Reference:  Tiwari, Jayaswal, Sinha (2021). Competitive Hub Location Problem: Model and Solution Approaches, Transportation Research-B, Under Review.


1. Fullcode_MA and Fullcode_SA implement the following three methods: Mixed Integer Second Order Conic Program (MISOCP), Cutting Plane Algorithm (CPA) and Lagrangian Relaxation within Cutting Plane Algorithm (LR-CPA).


2. These codes correspond to the multiple allocation (MA) version of the problem as described in Section 3.2 and the single allocation (SA) version of the problem as described in Section 3.1 of the manuscript.


3. The experiments are run on the following standard data-sets: CAB (10, 15, 20 and 25 nodes) and AP (50 nodes). Please refer to Section 6.1 of the manuscript. All the parameters used in obtaining the results as published in Section 6.2 are generated in the code.


4. The codes shared are in .cpp format and CPLEX 12.7.1 is used as the solver, where the default number of threads is 20. The code can be run using the following steps:
	
	i. Build the .exe file (using Microsoft Visual Studio). This generates "project name.exe", where "project name" is the name of the Microsoft Visual C++ project used to            build the .exe file.
	ii. Save a copy of the data file (e.g. Data_HS_CAB_10node.txt or AP_50nodes_data.txt) in the same folder as the .exe file
	iii. Open the Windows Command prompt. Change the directory to the location of the .exe file.
	iv. Run the following command on the command prompt: "project name.exe Data file name p alpha". Here p and alpha correspond to the parameter for which the code is to be             run.


5. To save the output to a file, run the following command on the command prompt: "project name.exe Data file name p alpha" > "Output file Name". 
