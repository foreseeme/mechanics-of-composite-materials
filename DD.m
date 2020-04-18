%求D矩阵时所需的函数
function F=DD(D,t,low,up)
n=up-low+1;
len=size(D,3);
F=zeros(3,3,len);
for j=1:3
    for i=1:3
        k=low;
        for kk=1:n
            F(i,j,kk)=D(i,j,kk)*((k*t)^3-((k-1)*t)^3)/3;
            kk=kk+1;
            k=k+1;
        end
        i=i+1;
    end
    j=j+1;
end