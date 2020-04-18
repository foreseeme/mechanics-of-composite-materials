%复合材料力学，强度问题计算，用于预测材料极限强度
%先输入参数
n=32
%Q=[141.9 3.06 0;3.06 8.66 0;0 0 5.0];单层刚度矩阵  HT3/5224材料
for ii=1:n
    Q(:,:,ii)=[141.9 3.06 0;3.06 8.66 0;0 0 5.0];   
end
S=zeros(3,3);%单层柔度矩阵
A=zeros(3,3);%层合板刚度矩阵
elt=zeros(3,1,n),exy=zeros(3,1,n),sgmlt=zeros(3,1,n),sgmxy=zeros(3,1,n);%sgm表示应力σ，lt表示主方向，xy表示非主方向
fai=[45 45 -45 -45 45 45 -45 -45 45 45 -45 -45 45 45 -45 -45 -45 -45 45 45 -45 -45 45 45 -45 -45 45 45 -45 -45 45 45];%主方向与xy方向的夹角
%m=cos(fai),n=sin(fai);
%T=[m*m n*n 2*m*n;n*n m*m -2*m*n;-m*n m*n m*m-n*n];%sgmlt=T*sgmxy
Qba=zeros(3,3,n);%非材料主方向单层刚度矩阵(k为层数)

Qba=QQ(Q,fai);
%S=inv(Q);暂时用不上S
%N，M为外载荷
N=[0;0;0];

M=[500;0;0];
exy0=zeros(1,3),kai=zeros(1,3);
t=4/32;
k=16;

%fai=[0,90,0,90,0,90,0,90,90,0,90,0,90,0,90,0];
aQba=AA(Qba,t,-k+1,k);
bQba=BB(Qba,t,-k+1,k);
dQba=DD(Qba,t,-k+1,k);

A=sumMAtrix(aQba,-k+1,k);
B=sumMAtrix(bQba,-k+1,k);
D=sumMAtrix(dQba,-k+1,k);

%inv(A)*b不如A\b，b*inv(A)不如b/A
B2=-A\B;
C2=B/A;
D2=D-B*inv(A)*B;
a=inv(A)-B2*inv(D2)*C2;
b=B2/inv(D2);
%c=-D2\C2;
c=b';
d=inv(D2);
exy0=a*N+b*M;
kai=c*N+d*M;
for kk=-k+1:k
   exy(:,:,kk+k)=exy0+(kk+kk-1)*t*kai/2; 
   sgmxy(:,:,kk+k)=Qba(:,:,kk+k)*exy(:,:,kk+k);
   m=cosd(fai(kk+k));
   n=sind(fai(kk+k));
   T=[m^2 n^2 2*m*n;n^2 m^2 -2*m*n;-m*n m*n m^2-n^2];
   sgmlt(:,:,kk+k)=T*sgmxy(:,:,kk+k);
   kk=kk+1; 
end













