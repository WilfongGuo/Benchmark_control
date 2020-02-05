function [ x,nd ] = control( dz,dN  )
% ********************����������Ϣ***********************
%���룺
%z���ڽӹ�ϵ
%�����
%nd��driver�ڵ����Ŀ
%driver��driver�ڵ�ļ���
% ********************simple example***********************
%z=[1 2;1 3;2 4];
%[ nd,driver ] = control( z );
% ********************�����ڽӾ������Ϣ***********************
%z=[z;z(:,2) z(:,1)];
%N=max(max(z));

uz=unique(dz);
[ind,z1]=ismember(dz(:,1),uz);
[ind,z2]=ismember(dz(:,2),uz);
z=[z1 z2];

N=max(max(z));

M=size(z);
T=M(1,1);
zz=zeros(N);
yy=[];
for i=1:T
zz(z(i,2),z(i,1))=1;
end
g=zz;
% *************�����������㷨��������������ƥ��***************
N=length(g);
Z=[];
Z(:,1)=z(:,1);Z(:,2)=z(:,2)+N;
NM=max(max(Z));

B=zeros(NM);[c,d]=size(z);

for i=1:c
    
    B(Z(i,1),Z(i,2))=1;
    B(Z(i,2),Z(i,1))=1;
    
end

k=sum(B);
isolited=find(k==0)';

B=sparse(B);
F = matching(B);
f=[];
f(:,1)=[1:N]';
f(:,2)=F(1:N,1);
PP=[];kk=1;
for i=1:length(f)
    if f(i,2)~=0
        PP(kk,1)=f(i,1);PP(kk,2)=f(i,2)-N;
        kk=kk+1;
    end
end
NV=[1:N]';
d=[setdiff(NV,PP(:,2))];
driver=unique(d);
nd=length(driver);
x=zeros(dN,1);
x(uz(driver),1)=1;

A_adjacent=zeros(N);
for i=1:size(z,1)
    
        A_adjacent(z(i,2),z(i,1))=1;
        A_adjacent(z(i,1),z(i,2))=1;
       
end

k=sum(A_adjacent)';
x(find(k==0))=0;
%**********************************************************
index=length(find(x(:,1)~=0));

end

