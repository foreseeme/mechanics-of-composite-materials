%求A矩阵时所需的函数
function F=AA(A,t,low,up)
n=up-low+1;
len=size(A,3);
F=zeros(3,3,len);
for j=1:3
    for i=1:3
        k=low;
        for kk=1:n
            F(i,j,kk)=A(i,j,kk)*(k*t-(k-1)*t);
            kk=kk+1;
            k=k+1;
        end
        i=i+1;
    end
    j=j+1;
end