clc
clear
%   $Id: main_Benchmark_control.m Created at 2020-02-05 21:25:22 $
%   by Weifeng Guo, Zhengzhou University, China
%   Copyright (c) 2019-2023 by School of Electrical Engineering, Zhengzhou University, 
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%Remaind: Please install gurobi before running our code
%Remaind: Please install gurobi before running our code
%Remaind: Please install gurobi before running our code

%**************Part 1:Input the information of samples and network information****
%**************sample information**************
%Example 1:Chu single cell data 


unzip('sc_time_course_data.zip') 
[data_main,~,~]=importdata('sc_time_course_data.csv');
[~,data_attribute,~]=importdata('SraRunTableTime.xlsx');


data=data_main.data;gene_list=data_main.textdata(2:end,1);
ref_data=data(:,find(data_attribute==1));
%Example 2:TCGA-BRCA cancer data (or other cancer datasets)

unzip('BRCA_normal.zip') 
unzip('BRCA_tumor.zip')
expression_tumor_fileName = 'BRCA_tumor.txt';
expression_normal_fileName = 'BRCA_normal.txt';


[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);tumor_data=tumor.data;
[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);normal_data=normal.data;
data=tumor_data;ref_data=normal_data;
%**************the network construction information****
%if Network_index=1,we use CSN; if Network_index=2,we use SSN
%if Network_index=3,we use SPCC; if Network_index=4,we use LIONESS
Network_method_index=1;
%Network_method_index=2;
%Network_method_index=3;
%Network_method_index=4;

%%**************Part 2:PDC outputs the predicted combinational drugs****
%Note that the input variable "ref_data" only is used by SSN,although it is a input variable in our function;

[ MMS,MDS,NCU,NCD ] = benchmark_control( data,ref_data,gene_list,Network_method_index );

%%**************Part 3:save the result****

save Benchmark_network_control_results
