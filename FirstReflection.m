FOV=90*pi/180; % field of view of the RX
xr=1; yr=1; zr=1;
xt=2; yt=4; zt=1;
n=1;  % radiation beam
Tf=1; % Filter Ttansmission
Tc=1; % Concentrator Ttansmission
Pt=1e-3; % optical power
Ar=1e-4; % collection area
dA= (5*5)*10^(-4);
S=3*10^(8);
ze=3;
ye1=0;
ye2=8;
xe3=0;
xe4=4;
i=0;
j=0;
k=0;
l=0;
m=0;
for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8
   i=i+1;

  %R=sqrt((xr-xt)^2+(yt-yr)^2+(zr-zt)^2);
  R1=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R2=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  R=R1+R2;

  theta1=acos(2/R1);
  phi1=theta1;
  theta2=acos(2/R2);
  phi2=theta2;

% Pr=((n+1)/(2*pi*R*R))*Tf*Tc*Pt*(cos(theta))^n*cos(phi);
  Pre=((n+1)*(n+1)/(4*(pi)^2*(R1)^2*(R2)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta1))^n*cos(phi1)*(cos(theta2))^n*cos(phi2);
    T=R/S;
    A(i)=Pre;
    B(i)=T;
    end

end


for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j=j+1;

    R3=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R4=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
    R5=R3+R4;

    theta3=acos(4/R3);
    phi3=acos(abs(zt-ze1)/R3);
    theta4=acos(1/R4);
    phi4=acos(abs(zr-ze1)/R4);

    Pre1=((n+1)*(n+1)/(4*(pi)^2*(R3)^2*(R4)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3))^n*cos(phi3)*(cos(theta4))^n*cos(phi4);
    T1=R5/S;
    C(j)=Pre1;
    D(j)=T1;
    end
end


for xe2=0.05:0.05:4    %Plane 2
    for ze2=1:0.05:3      %ye2=8
   k=k+1;

    R6=sqrt((xt-xe2)^2+(yt-ye2)^2+(zt-ze2)^2);
    R7=sqrt((xr-xe2)^2+(yr-ye2)^2+(zr-ze2)^2);
    R8=R6+R7;

    theta5=acos(6/R6);
    phi5=acos(abs(zt-ze2)/R6);
    theta6=acos(7/R7);
    phi6=acos(abs(zr-ze2)/R7);

    Pre2=((n+1)*(n+1)/(4*(pi)^2*(R6)^2*(R7)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta5))^n*cos(phi5)*(cos(theta6))^n*cos(phi6);
    T2=R8/S;
    E(k)=Pre2;
    F(k)=T2;
    end
end

for ye3=0.05:0.05:8    %Plane 3
    for ze3=1:0.05:3    %xe3=0
   l=l+1;

    R9=sqrt((xt-xe3)^2+(yt-ye3)^2+(zt-ze3)^2);
    R10=sqrt((xr-xe3)^2+(yr-ye3)^2+(zr-ze3)^2);
    R11=R9+R10;

    theta7=acos(2/R9);
    phi7=acos(abs(zt-ze3)/R9);
    theta8=acos(1/R10);
    phi8=acos(abs(zr-ze3)/R10);

    Pre3=((n+1)*(n+1)/(4*(pi)^2*(R9)^2*(R10)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta7))^n*cos(phi7)*(cos(theta8))^n*cos(phi8);
    T3=R11/S;
    G(l)=Pre3;
    H(l)=T3;
    end
end

for ye4=0.05:0.05:8    %Plane 4
    for ze4=1:0.05:3    %xe4=4
   m=m+1;

    R12=sqrt((xt-xe4)^2+(yt-ye4)^2+(zt-ze4)^2);
    R13=sqrt((xr-xe4)^2+(yr-ye4)^2+(zr-ze4)^2);
    R14=R12+R13;

    theta9=acos(2/R12);
    phi9=acos(abs(zt-ze4)/R12);
    theta10=acos(3/R13);
    phi10=acos(abs(zr-ze4)/R13);

    Pre4=((n+1)*(n+1)/(4*(pi)^2*(R12)^2*(R13)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta9))^n*cos(phi9)*(cos(theta10))^n*cos(phi10);
    T4=R14/S;
    W(m)=Pre4;
    X(m)=T4;
    end
end

%  b=B*10^9;
%  B=floor(b);
%
%   d=D*10^9;
%  D=floor(d);
%
%  f=F*10^9;
%  F=floor(f);
%
%
%  h=H*10^9;
%  H=floor(h);
%
%  x=X*10^9;
%  X=floor(x);
%
%
% Time=[B,D,F,H,X];
% Power=[A,C,E,G,W];




A(find(isnan(A)))=0;                  % Ceiling

  new_plos=zeros(1,80);

    for q=1:80

      bb=(floor((B*1e9)))-q;

       [rows,cols]=find(bb==0);

       if(cols>0)

       prq=sum(A(cols));

       new_plos(q)=new_plos(q)+prq;

       end

    end


C(find(isnan(C)))=0;                 %plane1

  new_plos1=zeros(1,80);

    for q=1:80

      dd=(floor((D*1e9)))-q;

       [rows,cols]=find(dd==0);

       if(cols>0)

       prq=sum(C(cols));

       new_plos1(q)=new_plos1(q)+prq;

       end

    end



E(find(isnan(E)))=0;                   %plane2

  new_plos2=zeros(1,80);

    for q=1:80

     ff=(floor((F*1e9)))-q;

       [rows,cols]=find(ff==0);

       if(cols>0)

       prq=sum(E(cols));

       new_plos2(q)=new_plos2(q)+prq;

       end

    end


G(find(isnan(G)))=0;                 %plane3

  new_plos3=zeros(1,80);

    for q=1:80

      hh=(floor((H*1e9)))-q;

       [rows,cols]=find(hh==0);

       if(cols>0)

       prq=sum(G(cols));

       new_plos3(q)=new_plos3(q)+prq;

       end

    end


W(find(isnan(W)))=0;                 %plane4

  new_plos4=zeros(1,80);

    for q=1:80

      xx=(floor((X*1e9)))-q;

       [rows,cols]=find(xx==0);

       if(cols>0)

       prq=sum(W(cols));

       new_plos4(q)=new_plos4(q)+prq;

       end

    end

    new_plos_total=new_plos+new_plos1+new_plos2+new_plos3+new_plos4;
