function [ new_x,index ] = Opti_weight_ncd( dz,NN  )
%non-lineaqr controllability of directed networks
%    min <xi>  and max sum<wixi>
%    wi-wj + nxi>=1
%    0<=wi<=n-1
%     xi=0,1
%   input:
%         z:network structure
%         scores:the scores of each node
%         lamda:balence parameters
%   Output:
%         lc_index:the number of driver nodes

%***********************solve the problem**************************************
%**********************MATLAB2014******************************

uz=unique(dz);
[ind,z1]=ismember(dz(:,1),uz);
[ind,z2]=ismember(dz(:,2),uz);
z=[z1 z2];
N=max(max(z));

lamda=0;
%N=max(max(z));
scores=ones(N,1);

N1=length(scores);
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,2*N1);

for i=1:N2
    
        A_adjacent(i,z(i,1))=1;
        A_adjacent(i,z(i,2))=-1;
        A_adjacent(i,N1+z(i,1))=N1;
       
end

f1=zeros(N1,1);
f2 = ones(N1,1)-lamda*scores;%
f=[f1;f2];

A=-A_adjacent;
b=-ones(size(A,1),1);

model.A = sparse(A);
model.obj = f';
model.rhs = b;
model.sense=[];
for i=1:size(A,1)
model.sense = [model.sense,'<'];;
end

model.modelsense = 'min';

lb=zeros(size(A,2),1);
ub1=(N1-1)*ones(N1,1);ub2=ones(N1,1);ub=[ub1;ub2];
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
x=x(N1+1:end,:);

N=max(max(z));
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N);
for i=1:size(z,1)
    
        A_adjacent(z(i,2),z(i,1))=1;
        %A_adjacent(z(i,1),z(i,2))=1;
       
end
k1=sum(A_adjacent)';
k2=sum(A_adjacent,2);
source=intersect(find(k1~=0),find(k2==0));
x(source,:)=1;


object=result.objval;
index=sum(x);

new_x=zeros(NN,1);
new_x(uz(find(x~=0)),1)=1;

end

