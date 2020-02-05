# Instruction of Benchmark_control package
%Remaind: Please install gurobi before running our code (http://www.gurobi.com/) 

%Remaind: Please install gurobi before running our code (http://www.gurobi.com/) 

%Remaind: Please install gurobi before running our code (http://www.gurobi.com/) 

This package includes Matlab scripts and several datasets for demo of network control approach: 

(a)	main_Benchmark_control.m is a Matlab function for the routine of experimental analysis. 

(b) benchmark_control.m is the main script to call Benchmark_control 

(c) The input datasets include: % data:the tumor expression data % gene_list:the gene list name data
% ref_data:the reference data used in SSN 

% index:denotes we use which network construction method 

%if index=1,we use CSN 

%if index=2,we use SSN 

%if index=3,we use SPCC 

%if index=2,we use LIONESS

The output datasets include: The sample-specific driver profiles (matrix) by using MMS,MDS,NCU,NCD； For “MMS or MDS,NCU,NCD”, the column is the samples and the rows is the genes. The value “1” denoted that the gene is driver genes; 

(d) As a demo, users can directly run main_Benchmark_control.m in Matlab. We choose the single cell time cource data and BRCA cancer data as a test case in our demo. This package has been tested in different computer environments as: Window 7 or above; Matlab 2014 or above.

(e) When users analyzed yourself new data, please: (1) Prepare input datasets as introduced in (d). (2) Clear the previous results. (3) Set parameters in benchmark_control.m as introduced in (b). (4) Run main_Benchmark_control.m. (5) Suggest that the users add all fille in our folders to your folder.

% $Id: main_Benchmark_control.m Created by Weifeng Guo, Zhengzhou University, China at 2020-02-05 21:25:22 $ % $Copyright (c) 2019-2022 by School of Electrical Engineering, Zhengzhou University, Zhengzhou$;

% $If any problem, pleasse contact shaonianweifeng@126.com for help. $
