%求B矩阵时所需的函数
function F=BB(B,t,low,up)
n=up-low+1;
len=size(B,3);
F=zeros(3,3,len);
for j=1:3
    for i=1:3
        k=low;
        for kk=1:n
            F(i,j,kk)=B(i,j,kk)*((k*t)^2-((k-1)*t)^2)/2;
            kk=kk+1;
            k=k+1;
        end
        i=i+1;
    end
    j=j+1;
end