function [ MMS,MDS,NCU,NCD ] = benchmark_control( data,ref_data,gene_list,Network_method_index )
%we output the sample-specific driver profiles by using different control
%methods
%   Input:
%         data:the tumor expression data
%         gene_list:the gene list name data  
%         ref_data:the reference data used in SSN
%         index=1 or index=2:denotes we use which network construction method
%**************the network information****
%if index=1,we use CSN
%if index=2,we use SSN
%if Network_index=3,we use SPCC; 
%if Network_index=2,we use LIONESS
%  Output:
%         The sample-specific driver profiles of MMS,MDS,NCU,NCD;
%         The column is the samples and the rows is the genes. The value “1” denoted that the gene is driver genes; 
%************************part1:LOAD sample data and network data************************
%********************obtain the paired expression data******************
index=Network_method_index;
load('GIN_network_information.mat')% Hou JP, Ma J. DawnRank: discovering personalized driver genes in cancer. Genome Medicine. 2014;6(7):56. doi: 10.1186/S14073-014-0056-8.
%load('sPPI_network_information.mat') % Vinayagam A, Stelzl U, Foulle R, Plassmann S, Zenkner M, Timm J, et al. A Directed Protein Interaction Network for Investigating Intracellular Signal Transduction. Science Signaling. 2011;4(189):rS9-rs. doi: 10.1126/scisignal.2001699.
%[~,PPI,~]=xlsread('network_FIsInGene_041709.xlsx'); % Bertrand D, Chng KR, Sherbaf FG, Kiesel A, Chia BK, Sia YY, et al. Patient-specific driver gene prediction and risk assessment through integrated network analysis of cancer omics profiles. Nucleic acids research 2015;43:e44-e
[x1,y1]=ismember(edge0(:,1),gene_list);
[x2,y2]=ismember(edge0(:,2),gene_list);

y=y1.*y2;
z=[y1 y2];
z(find(y==0),:)=[];
N1=length(gene_list);
[N2,~]=size(z);

Net=zeros(N1,N1);

for i=1:N2
    
         Net(z(i,2),z(i,1))=1;  %undirected gene-gene interaction network
         %Net(z(i,1),z(i,2))=1;    
end

%ndm = csndm(a,0.01,0.1,0);

if index==1
    
 %***************CSN*********************** 
  tic
  
for i=1:size(data,2)
    
   
    i
    csn = csnet(data,i,0.01,0.1,0);
    cand=full(csn{1,i});
    C=cand.*Net;
    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

    
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
NCD(:,i)=NCD_x1;
NCU(:,i)=NCU_x1;
    
    
 
    
end

toc

end


if index==2
    
    tic
%***************SSN*********************** 
 
for i=1:size(data,2)
    i
    sample=data(:,i);
    [ index_R,p ] =SSN(sample,ref);
    p(isnan(p))=1;
    p(p>=0.01)=0;
    p(p~=0)=1;
    
    cand=p.*Net;
    
    C=cand;
    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

    
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
NCD(:,i)=NCD_x1;
NCU(:,i)=NCU_x1;


    clear index_R p cand

end


toc

end


if index==3
    
%***************SPCC***********************   
    tic
  for i=1:size(data,2)  
    i
    csn = spcc_method(data,i);
    cand=csn;

    C=cand.*Net;
    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

    
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
NCD(:,i)=NCD_x1;
NCU(:,i)=NCU_x1;
    
clear C
    
  end

    toc
    
end

if index==4
    
%***************LIONESS***********************   
   tic 
   for i=1:size(data,2)
   i
    p=lioness_method(data,i);
  
    C=p.*Net; 

    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

   
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
NCD(:,i)=NCD_x1;
NCU(:,i)=NCU_x1;
    
clear C
    
   
    
   end
   
   
toc
   
end




end

