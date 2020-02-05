function [ new_x,index ] = Opti_MDS( dz,NN )
%MDS controllability of undirected networks
%   input:
%         z:network structure
%   Output:
%         lc_index:the number of driver nodes

%***********************solve the problem**************************************
%**********************MATLAB2014******************************
%N=max(max(z));
%calculate the adjacency matrix of bipartite graph


uz=unique(dz);
[ind,z1]=ismember(dz(:,1),uz);
[ind,z2]=ismember(dz(:,2),uz);
z=[z1 z2];
N=max(max(z));



A_adjacent=zeros(N);
for i=1:size(z,1)
    
        A_adjacent(z(i,2),z(i,1))=1;
        A_adjacent(z(i,1),z(i,2))=1;
       
end

k=sum(A_adjacent)';
A=-(A_adjacent+eye(N));
b=-ones(size(A,1),1);
f = ones(size(A,2),1);%

model.A = sparse(A);
model.obj = f';
model.rhs = b;
model.sense=[];
for i=1:size(A,1)
model.sense = [model.sense,'<'];;
end

model.modelsense = 'min';


lb=zeros(size(A,2),1);
ub=ones(size(A,2),1);

 model.lb = lb;
 model.ub = ub;

n = size(model.A, 2);


model.vtype = repmat('C', n, 1);

intcon = [1:size(A,2)];
model.vtype(intcon) = 'I';
%model.varnames = names;

%gurobi_write(model, 'mip1.lp');

params.outputflag = 0;

result = gurobi(model, params);

%disp(result);
x=result.x;
x(find(k==0))=0;


object=result.objval;
index=sum(x);

new_x=zeros(NN,1);
new_x(uz(find(x~=0)),1)=1;




end

