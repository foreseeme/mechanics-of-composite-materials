function F=sumMAtrix(A,low,up)
%size(A,3)可以得到A的第三维的长度
F=zeros(3,3);
n=up-low+1;
for i=1:n
    F=A(:,:,i)+F;   
end