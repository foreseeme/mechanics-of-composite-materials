%º∆À„F
function F=QQ(Q,fai)
nn=length(fai);

F=zeros(3,3,nn);
for k=1:nn 
    m=cosd(fai(k));
    n=sind(fai(k));
    F(1,1,k)=m^4*Q(1,1,k)+n^4*Q(2,2,k)+2*m*m*n*n*Q(1,2,k)+4*m*m*n*n*Q(3,3,k);
    F(2,2,k)=n^4*Q(1,1,k)+m^4*Q(2,2,k)+2*m*m*n*n*Q(1,2,k)+4*m*m*n*n*Q(3,3,k);
    F(1,2,k)=m*m*n*n*Q(1,1,k)+m*m*n*n*Q(2,2,k)+(m^4+n^4)*Q(1,2,k)+(-4*m*m*n*n)*Q(3,3,k);
    F(3,3,k)=m*m*n*n*Q(1,1,k)+m*m*n*n*Q(2,2,k)+(-2*m*m*n*n)*Q(1,2,k)+(m*m-n*n)^2*Q(3,3,k);
    F(1,3,k)=m*m*m*n*Q(1,1,k)-m*n*n*n*Q(2,2,k)+(m*n*n*n-m*m*m*n)*Q(1,2,k)+2*(m*n*n*n-m*m*m*n)*Q(3,3,k);
    F(2,3,k)=m*n*n*n*Q(1,1,k)-m*m*m*n*Q(2,2,k)+(m*m*m*n-m*n*n*n)*Q(1,2,k)+2*(m*m*m*n-m*n*n*n)*Q(3,3,k);
    F(2,1,k)=F(1,2,k);
    F(3,1,k)=F(1,3,k);
    F(3,2,k)=F(2,3,k);
    k=k+1;
end
                
