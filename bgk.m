%The aim of this function is to estimate ||bgk||^2 vector
function [bgk2]=bgk(b,K,group_size)
bgk2=zeros(K,1);
for i=1:K
   bgk2(i,1)=b(group_size==i)'*b(group_size==i); 
end
end