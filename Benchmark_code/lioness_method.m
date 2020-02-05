function [CSN] = lioness_method(data,i)
%construct sample specific network


PCC=corrcoef(data');
v0=data;

    v0(:,i)=[];
    PCC1=corrcoef(v0');
    x=size(data,2)*(PCC-PCC1)+PCC1;
    x(isnan(x))=0;
    
    xx=triu(x);
    a_x=abs(xx(xx~=0));
    threshold=mean(a_x)+2*std(a_x);
    
    x(abs(x)<threshold)=0;
    x(x~=0)=1;
  
    CSN=x;




end

