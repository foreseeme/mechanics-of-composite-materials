%蔡吴准则
clear
%层合板层数
n=16;
%Q=[141.9 3.06 0;3.06 8.66 0;0 0 5.0];HT3/5224材料的Q矩阵
Q=zeros(3,3,n);
for ii=1:n
    Q(:,:,ii)=[141.9 3.06 0;3.06 8.66 0;0 0 5.0];  %每层的Q矩阵 
end
Xt=1400;Xc=1100;Yt=50;Yc=180;S=99;
F11=1/(Xt*Xc);F22= 1/(Yt*Yc);F33=1/(S^2);F1=1/Xt-1/Xc;F2=1/Yt-1/Yc;
%F12=-0.5*sqrt(F11*F22);
F12=-2.25*10^(-6);
X=0;Y=0;S=99;%蔡-希尔准则参数?
elt=zeros(3,1,n);exy=zeros(3,1,n);sgmlt=zeros(3,1,n);sgmxy=zeros(3,1,n);%sgmlt表示主方向应力sgmx表示正方向应力
fai=[0 90 0 90 0 90 0 90 90 0 90 0 90 0 90 0];%各层角度?
Qba=zeros(3,3,n);
N=[0;0;0];
Mx=100;
M=[Mx;0;0];%N的单位是N/mm，M的单位是N*mm
exy0=zeros(3,1);kai=zeros(3,1);
t=0.125;
k=8;
R=zeros(n,1);
WW=0;%极限强度
xiabiao=0;
a=1;
while (a==1)
   Qba=QQ(Q,fai);
   aQba=AA(Qba,t,-k+1,k);
   bQba=BB(Qba,t,-k+1,k);
   dQba=DD(Qba,t,-k+1,k);
   
   A=sumMAtrix(aQba,-k+1,k);
   B=sumMAtrix(bQba,-k+1,k);
   D=sumMAtrix(dQba,-k+1,k);
  
    %inv(A)*b可用A\b代替，b*inv(A)可用b/A代替
   B2=-A\B;
   C2=B/A;
   D2=D-B*inv(A)*B;
   aa=inv(A)-B2*inv(D2)*C2;
   bb=B2*inv(D2);
   cc=-D2\C2;
  % cc=bb';
   dd=inv(D2);
   exy0=(aa*N+bb*M);
   kai=(cc*N+dd*M);

   for kk=-k+1:1:k       
       exy(:,:,kk+k)=exy0+(kk+kk-1)*t*kai/2;
       sgmxy(:,:,kk+k)=Qba(:,:,kk+k)*exy(:,:,kk+k);       
       m=cosd(fai(kk+k));
       n=sind(fai(kk+k));
       T=[m^2 n^2 2*m*n;n^2 m^2 -2*m*n;-m*n m*n m^2-n^2];
       sgmlt(:,:,kk+k)=T*sgmxy(:,:,kk+k);
       %Two表示二次项，One表示一次项
       Two=F11*sgmlt(1,1,kk+k)^2+F22*sgmlt(2,1,kk+k)^2+F33*sgmlt(3,1,kk+k)^2+2*F12*sgmlt(1,1,kk+k)*F22*sgmlt(2,1,kk+k);
       One=F1*sgmlt(1,1,kk+k)+F2*sgmlt(2,1,kk+k);
       R(kk+k)=(-One+sqrt(One^2+4*Two))/(2*Two);
   end
   [m,index]=min(R);
   El=140;Et=8.6*0.3;Glt=5.0*0.3;niult=0.3*0.35;niutl=niult*Et/El;
   Q(1,1,index)=El/(1-niult*niutl);
   Q(1,2,index)=niutl*El/(1-niult*niutl);Q(2,1,index)=Q(1,2,index);
   Q(2,2,index)=Et/(1-niult*niutl);
   Q(3,3,index)=Glt;
   Q(1,3,index)=0; Q(3,1,index)=0; Q(2,3,index)=0; Q(3,2,index)=0;
   WWW=R(index)*Mx;
   if WWW>WW
       WW=WWW;
       xiabiao=index;
  
   else
       a=0;       
   end
   
end