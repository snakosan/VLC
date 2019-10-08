
clear all;
FOV=90*pi/180; % field of view of the RX
xr=1; yr=7; zr=1;
xt=2; yt=4; zt=1;
n=1;  % radiation beam
Tf=1; % Filter Ttansmission
Tc=1; % Concentrator Ttansmission
Pt=1; % optical power
Ar=1e-4; % collection area
P1=0.8;
P2=0.8;
dA=(5*5)*10^(-4);
dA2=(20*20)*10^(-4);
S=3*10^(8);
ze=3;
ye1=0;
ye2=8;
xe3=0;
xe4=4;
i=0;
i_1=0;
i_2=0;
i_3=0;
i_4=0;
i_5=0;
i_6=0;   
i_7=0;
j=0;
j_1=0;
j_2=0;
j_3=0;
j_4=0;
j_5=0;
j_6=0;
j_7=0;
j_8=0;
k=0;
k_1=0;
k_2=0;
k_3=0;
k_4=0;
k_5=0;
k_6=0;
k_7=0;
k_8=0;
l=0;
l_1=0;
l_2=0;
l_3=0;
l_4=0;
l_5=0;
l_6=0;
l_7=0;
l_8=0;
m=0;
m_1=0;
m_2=0;
m_3=0;
m_4=0;
m_5=0;
m_6=0;
m_7=0;
m_8=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % FIRST REFLECTION
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8
   i=i+1;
       
  %R=sqrt((xr-xt)^2+(yt-yr)^2+(zr-zt)^2);
  R1=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R2=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  R=R1+R2;

  theta1=acos(abs(ze-zt)/R1);         % ze-zt=3-1=2
  phi1=theta1;
  theta2=acos(abs(ze-zt)/R2);         % ze-zr=3-1=2
  phi2=theta2;

% Pr=((n+1)/(2*pi*R*R))*Tf*Tc*Pt*(cos(theta))^n*cos(phi);
  Pre=((n+1)*(n+1)/(4*(pi)^2*(R1)^2*(R2)^2))*Tf*Tc*Pt*dA*P1*Ar*(cos(theta1))^n*cos(phi1)*(cos(theta2))^n*cos(phi2);
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
    
    theta3=acos(abs(yt-ye1)/R3);              % yt-ye1=4-0=4
    phi3=acos(abs(zt-ze1)/R3);
    theta4=acos(abs(yr-ye1)/R4);               %yr-ye1=1-0=1 
    phi4=acos(abs(zr-ze1)/R4);
    
    Pre1=((n+1)*(n+1)/(4*(pi)^2*(R3)^2*(R4)^2))*Tf*Tc*Pt*dA*P1*Ar*(cos(theta3))^n*cos(phi3)*(cos(theta4))^n*cos(phi4);
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
    
    theta5=acos(abs(yt-ye2)/R6);            % yt-ye2=4-8=4
    phi5=acos(abs(zt-ze2)/R6);
    theta6=acos(abs(yr-ye1)/R7);            % yr-ye1=1-8=7
    phi6=acos(abs(zr-ze2)/R7);
    
    Pre2=((n+1)*(n+1)/(4*(pi)^2*(R6)^2*(R7)^2))*Tf*Tc*Pt*dA*P1*Ar*(cos(theta5))^n*cos(phi5)*(cos(theta6))^n*cos(phi6);
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
    
    theta7=acos(abs(xt-xe3)/R9);      %xt-xe3=2-0=2
    phi7=acos(abs(zt-ze3)/R9);
    theta8=acos(abs(xt-xe3)/R10);      %xr-xt3=1-0=1
    phi8=acos(abs(zr-ze3)/R10);
    
    Pre3=((n+1)*(n+1)/(4*(pi)^2*(R9)^2*(R10)^2))*Tf*Tc*Pt*dA*P1*Ar*(cos(theta7))^n*cos(phi7)*(cos(theta8))^n*cos(phi8);
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
    
    theta9=acos(abs(xt-xe4)/R12);       % xe4-xt=4-2=2
    phi9=acos(abs(zt-ze4)/R12);
    theta10=acos(abs(xr-xe4)/R13);      % xe4-xr=4-1=3
    phi10=acos(abs(zr-ze4)/R13);
    
    Pre4=((n+1)*(n+1)/(4*(pi)^2*(R12)^2*(R13)^2))*Tf*Tc*Pt*dA*P1*Ar*(cos(theta9))^n*cos(phi9)*(cos(theta10))^n*cos(phi10);
    T4=R14/S;
    W(m)=Pre4;
    X(m)=T4;
    end
end




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

   
   FirstReflection_plos=new_plos+new_plos1+new_plos2+new_plos3+new_plos4;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SECOND REFLECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ( Ceiling with all Planes )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 

for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3                                                                                               
   i=i+1;
     
   for xe1=0.2:0.2:4     %Plane 1
    for ze1=1:0.2:3        %ye1=0
   i_1=i_1+1;
   
   R1=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R3=sqrt((xe-xe1)^2+(ye-ye1)^2+(ze-ze1)^2);
   R4=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R=R1+R3+R4;
  
  %angel for Second Ref
  theta11=acos(2/R1);
  theta12=acos(2/R1);            
  theta13=acos(abs(ye-ye1)/R3);             
  theta14=acos(abs(ze-ze1)/R3);   
  theta15=acos(abs(yr-ye1)/R4);             %yr-ye1=1-0=1
  theta16=acos(abs(zr-ze1)/R4); 
  
  Pres=(((n+1)*(n+1))^2/(8*(pi)^3*(R1)^2*(R3)^2*(R4)^2))*Tf*Tc*Pt*dA*P1*dA2*P2*Ar*(cos(theta11))^n*(cos(theta12))^n*cos(theta13)*(cos(theta14))^n*cos(theta15)*cos(theta16);
   T=R/S;
   B(i_1)=Pres;
   C(i_1)=T;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3
   i_2=i_2+1;
  
 
   for xe2=0.2:0.2:4     %Plane 2
    for ze2=1:0.2:3        %ye2=8
      i_3=i_3+1;
   
   R11=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R13=sqrt((xe-xe2)^2+(ye-ye2)^2+(ze-ze2)^2);
   R14=sqrt((xr-xe2)^2+(yr-ye2)^2+(zr-ze2)^2);
   R15=R11+R13+R14;
  
   %angel for Second Ref
  theta17=acos(2/R11);
  theta18=acos(2/R11);             
  theta19=acos(abs(ye-ye2)/R13);             
  theta20=acos(abs(ze-ze2)/R13);   
  theta21=acos(abs(yr-ye2)/R14);        %yr-ye2=1-8=7      
  theta22=acos(abs(zr-ze2)/R14);              
  
  Pres1=(((n+1)*(n+1))^2/(8*(pi)^3*(R11)^2*(R13)^2*(R14)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta17))^n*(cos(theta18))^n*cos(theta19)*(cos(theta20))^n*cos(theta21)*cos(theta22);
   T1=R15/S;
   B_1(i_3)=Pres1; 
   C_1(i_3)=T1;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3
   i_4=i_4+1;
  
   for ye3=0.2:0.2:8    %Plane 3
     for ze3=1:0.2:3    %xe3=0
           i_5=i_5+1;

     R16=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
     R18=sqrt((xe-xe3)^2+(ye-ye3)^2+(ze-ze3)^2);
     R19=sqrt((xr-xe3)^2+(yr-ye3)^2+(zr-ze3)^2);
     R20=R16+R18+R19;

        %angel for Second Ref
      theta23=acos(2/R16);
      theta24=acos(2/R16);             
      theta25=acos(abs(xe-xe3)/R18);            
      theta26=acos(abs(ze-ze3)/R18);
      theta27=acos(abs(xr-xe3)/R19);             %xr-xe3=1-0=1
      theta28=acos(abs(zr-ze3)/R19);           

    Pres2=(((n+1)*(n+1))^2/(8*(pi)^3*(R16)^2*(R18)^2*(R19)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta23))^n*(cos(theta24))^n*cos(theta25)*(cos(theta26))^n*cos(theta27)*cos(theta28);
         T2=R20/S;
         B_2(i_5)=Pres2; 
         C_2(i_5)=T2;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3                                                                                               
   i_6=i_6+1;
   
   for ye4=0.2:0.2:8    %Plane 4
     for ze4=1:0.2:3    %xe4=4
   i_7=i_7+1;
   
   R21=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R23=sqrt((xe-xe4)^2+(ye-ye4)^2+(ze-ze4)^2);
   R24=sqrt((xr-xe4)^2+(yr-ye4)^2+(zr-ze4)^2);
   R25=R21+R23+R24;
  
    %angel for Second Ref
  theta29=acos(2/R21);
  theta30=acos(2/R21);             
  theta31=acos(abs(xe-xe4)/R23);       
  theta32=acos(abs(ze-ze4)/R23);
  theta33=acos(abs(xr-xe4)/R24);             %xr-xe4=1-4=3
  theta34=acos(abs(zr-ze4)/R24);             
  
  Pres3=(((n+1)*(n+1))^2/(8*(pi)^3*(R21)^2*(R23)^2*(R24)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta29))^n*(cos(theta30))^n*cos(theta31)*(cos(theta32))^n*cos(theta33)*cos(theta34);
    T3=R25/S;
   
   B_3(i_7)=Pres3; 
   C_3(i_7)=T3;
    end
       end
    end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ( plane 1 with all )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j_1=j_1+1;
    
   for xe=0.2:0.2:4      %ceiling
   for ye=0.2:0.2:8      %ze=3                                                                                               
   j_2=j_2+1;
   
   R26=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
   R28=sqrt((xe-xe1)^2+(ye-ye1)^2+(ze-ze1)^2);
   R29=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R30=R26+R28+R29;
    
     %angel for Second Ref
  theta37=acos(abs(zt-ze1)/R26);
  theta38=acos(abs(yt-ye1)/R26);             %yt-ye1=4-0=4
  theta39=acos(abs(ze-ze1)/R28);
  theta40=acos(abs(ye-ye1)/R28);    
  theta41=acos(2/R29);            
  theta42=acos(2/R29);     
     
  Pres4=(((n+1)*(n+1))^2/(8*(pi)^3*(R26)^2*(R28)^2*(R29)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta37))^n*(cos(theta38))^n*cos(theta39)*(cos(theta40))^n*cos(theta41)*cos(theta42);
    T4=R30/S;
    F_1(j_2)=Pres4;
    G_1(j_2)=T4;
    
    end
       end
    end
   
end

for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j_3=j_3+1;
    
   for xe2=0.2:0.2:4     %Plane 2
    for ze2=1:0.2:3        %ye2=8                                                                                             
   j_4=j_4+1;
   
   R31=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
   R33=sqrt((xe1-xe2)^2+(ye1-ye2)^2+(ze1-ze2)^2);
   R34=sqrt((xr-xe2)^2+(yr-ye2)^2+(zr-ze2)^2);
   R35=R31+R33+R34;
    
   %angel for Second Ref
  theta45=acos(abs(zt-ze1)/R31);
  theta46=acos(abs(yt-ye1)/R31);        %ye1-yt=4-0=4
  theta47=acos(abs(ye2-ye1)/R33);      %%OR %ye2-ye1=8-0=8
  theta48=acos(abs(ze2-ze)/R33);
  theta49=acos(abs(zr-ze2)/R34);           
  theta50=acos(abs(yr-ye2)/R34);         %ye2-yr=8-1=7    
  
     
  Pres5=(((n+1)*(n+1))^2/(8*(pi)^3*(R31)^2*(R33)^2*(R34)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta45))^n*(cos(theta46))^n*cos(theta47)*(cos(theta48))^n*cos(theta49)*cos(theta50);
    T5=R35/S;
    F_2(j_4)=Pres5;
    G_2(j_4)=T5;
    
    end
       end
    end
   
end

for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j_5=j_5+1;
   
     for ye3=0.2:0.2:8    %Plane 3
     for ze3=1:0.2:3    %xe3=0
           j_6=j_6+1;

      R36=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
      R38=sqrt((xe1-xe3)^2+(ye1-ye3)^2+(ze1-ze3)^2);
      R39=sqrt((xr-xe3)^2+(yr-ye3)^2+(zr-ze3)^2);
      R40=R36+R38+R39;

      %angel for Second Ref
      theta53=acos(abs(zt-ze1)/R36);
      theta54=acos(abs(yt-ye1)/R36);            %yt-ye1=4-0=4
      theta55=acos(abs(xe1-xe3)/R38);        
      theta56=acos(abs(ze1-ze3)/R38);
      theta57=acos(abs(zr-ze3)/R39);            
      theta58=acos(abs(yr-ye3)/R39);         

     Pres6=(((n+1)*(n+1))^2/(8*(pi)^3*(R36)^2*(R38)^2*(R39)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta53))^n*(cos(theta54))^n*cos(theta55)*(cos(theta56))^n*cos(theta57)*cos(theta58);
        T6=R40/S;
        F_3(j_6)=Pres6; 
        G_3(j_6)=T6;
    end
       end
    end
   
end

   for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j_7=j_7+1;
   
       for ye4=0.2:0.2:8    %Plane 4
          for ze4=1:0.2:3    %xe4=4
             j_8=j_8+1;
   
    R41=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R43=sqrt((xe1-xe4)^2+(ye1-ye4)^2+(ze1-ze4)^2);
    R44=sqrt((xr-xe4)^2+(yr-ye4)^2+(zr-ze4)^2);
    R45=R41+R43+R44;
   
       %angel for Second Ref
    theta61=acos(abs(zt-ze1)/R41);
    theta62=acos(abs(yt-ye1)/R41);             %yt-ye1=4-0=4
    theta63=acos(abs(xe1-xe4)/R43);         
    theta64=acos(abs(ze1-ze4)/R43);
    theta65=acos(abs(zr-ze4)/R44);             
    theta66=acos(abs(yr-ye4)/R44);             
  
  Pres7=(((n+1)*(n+1))^2/(8*(pi)^3*(R41)^2*(R43)^2*(R44)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta61))^n*(cos(theta62))^n*cos(theta63)*(cos(theta64))^n*cos(theta65)*cos(theta66);
    T7=R45/S;
   F_4(j_8)=Pres7; 
   G_4(j_8)=T7;
    end
       end
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ( plane 2 with all )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for xe2=0.05:0.05:4        %Plane 2
    for ze2=1:0.05:3        %ye2=8
      k_1=k_1+1;
   
    for xe=0.2:0.2:4       %ceiling
   for ye=0.2:0.2:8        %ze=3                                                                                               
     k_2=k_2+1;
     
   R46=sqrt((xt-xe2)^2+(yt-ye2)^2+(zt-ze2)^2);
   R48=sqrt((xe2-xe)^2+(ye2-ye)^2+(ze2-ze)^2);
   R49=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R50=R46+R48+R49;
    
     %angel for Second Ref
  theta69=acos(abs(yt-ye2)/R46);      %ye2-yt=8-4=6
  theta70=acos(abs(zt-ze2)/R46);
  theta71=acos(abs(ze-ze2)/R48);
  theta72=acos(abs(ye-ye2)/R48);       %ye2-ye 
  theta73=acos(2/R49);       
  theta74=acos(2/R49);     
     
  Pres8=(((n+1)*(n+1))^2/(8*(pi)^3*(R46)^2*(R48)^2*(R49)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta69))^n*(cos(theta70))^n*cos(theta71)*(cos(theta72))^n*cos(theta73)*cos(theta74);
    T8=R50/S;
    I_1(k_2)=Pres8;
    H_1(k_2)=T8;
    
    end
       end
    end
   
end
    
for xe2=0.05:0.05:4          %Plane 2
    for ze2=1:0.05:3           %ye2=8
      k_3=k_3+1;
      
      for xe1=0.2:0.2:4       %Plane 1
    for ze1=1:0.2:3              %ye1=0
      k_4=k_4+1;
   
   R51=sqrt((xt-xe2)^2+(yt-ye2)^2+(zt-ze2)^2);
   R53=sqrt((xe2-xe1)^2+(ye2-ye1)^2+(ze2-ze1)^2);
   R54=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R55=R51+R53+R54;
  
  %angel for Second Ref
  theta77=acos(abs(yt-ye2)/R51);            %ye2-yt=8-4=4
  theta78=acos(abs(zt-ze2)/R51);
  theta79=acos(abs(ye1-ye2)/R53);           %%OR %ye1-ye2=8-0=8 // ze1-ze2
  theta80=acos(abs(ze-ze1)/R53);
  theta81=acos(abs(yr-ye1)/R54);             %yr-ye1=1-0=1
  theta82=acos(abs(zr-ze1)/R54);       
  
  Pres9=(((n+1)*(n+1))^2/(8*(pi)^3*(R51)^2*(R53)^2*(R54)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta77))^n*(cos(theta78))^n*cos(theta79)*(cos(theta80))^n*cos(theta81)*cos(theta82);
   T9=R55/S;
   I_2(k_4)=Pres9;
   H_2(k_4)=T9;
    end
       end
    end  
end

for xe2=0.05:0.05:4          %Plane 2
 for ze2=1:0.05:3             %ye2=8
  k_5=k_5+1;

    for ye3=0.2:0.2:8       %Plane 3
     for ze3=1:0.2:3          %xe3=0
      k_6=j_6+1;

     R56=sqrt((xt-xe2)^2+(yt-ye2)^2+(zt-ze2)^2);
     R58=sqrt((xe2-xe3)^2+(ye2-ye3)^2+(ze2-ze3)^2);
     R59=sqrt((xr-xe3)^2+(yr-ye3)^2+(zr-ze3)^2);
     R60=R56+R58+R59;
     
       %angel for Second Ref
      theta85=acos(abs(yt-ye2)/R56);           %ye2-yt=8-4=4
      theta86=acos(abs(zt-ze2)/R56);
      theta87=acos(abs(xe2-xe3)/R58);      
      theta88=acos(abs(ze2-ze3)/R58);
      theta89=acos(abs(zr-ze3)/R59);   
      theta90=acos(abs(yr-ye3)/R59);            
      
      Pres10=(((n+1)*(n+1))^2/(8*(pi)^3*(R56)^2*(R58)^2*(R59)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta85))^n*(cos(theta86))^n*cos(theta87)*(cos(theta88))^n*cos(theta89)*cos(theta90);
  
          T10=R60/S;
          I_3(k_6)=Pres10; 
          H_3(k_6)=T10;
 
     end
    end
   end
end
 

for xe2=0.05:0.05:4         %Plane 2
 for ze2=1:0.05:3            %ye2=8
  k_7=k_7+1;

       for ye4=0.2:0.2:8        %Plane 4
        for ze4=1:0.2:3            %xe4=4
          k_8=k_8+1;
   
    R61=sqrt((xt-xe2)^2+(yt-ye2)^2+(zt-ze2)^2);
    R63=sqrt((xe2-xe4)^2+(ye2-ye4)^2+(ze2-ze4)^2);
    R64=sqrt((xr-xe4)^2+(yr-ye4)^2+(zr-ze4)^2);
    R65=R61+R63+R64;
  
      %angel for Second Ref
    theta93=acos(abs(yt-ye2)/R61);         %ye2-yt=8-4=4
    theta94=acos(abs(zt-ze2)/R61);
    theta95=acos(abs(xe2-xe4)/R63);             
    theta96=acos(abs(ze2-ze4)/R63);
    theta97=acos(abs(zr-ze4)/R64);    
    theta98=acos(abs(yr-ye4)/R64);              
  
    Pres11=(((n+1)*(n+1))^2/(8*(pi)^3*(R61)^2*(R63)^2*(R64)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta93))^n*(cos(theta94))^n*cos(theta95)*(cos(theta96))^n*cos(theta97)*cos(theta98);
     
      T11=R65/S; 
      I_4(k_8)=Pres11; 
      H_4(k_8)=T11;
   
     end
       end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         % ( plane 3 with all )

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 
for ye3=0.05:0.05:8       %Plane 3
    for ze3=1:0.05:3       %xe3=0
   l_1=l_1+1;
   
     for xe=0.2:0.2:4        %ceiling
      for ye=0.2:0.2:8         %ze=3                                                                                               
        l_2=l_2+1;
   
     R66=sqrt((xt-xe3)^2+(yt-ye3)^2+(zt-ze3)^2);
     R68=sqrt((xe-xe3)^2+(ye-ye3)^2+(ze-ze3)^2);
     R69=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
        R70=R66+R68+R69;
    
        %angel for Second Ref
     theta3_3=acos(abs(xt-xe3)/R66);             %xt-xe3=2-0=2  
     theta3_4=acos(abs(zt-ze3)/R66);
     theta3_5=acos(abs(ze-ze3)/R68);
     theta3_6=acos(abs(xe-xe3)/R68);             % OR ye3-ye  
     theta3_7=acos(2/R69);             
     theta3_8=acos(2/R69);     
     
      Pres12=(((n+1)*(n+1))^2/(8*(pi)^3*(R66)^2*(R68)^2*(R69)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta3_3))^n*(cos(theta3_4))^n*cos(theta3_5)*(cos(theta3_6))^n*cos(theta3_7)*cos(theta3_8);
   
       T12=R70/S;
       K_1(l_2)=Pres12;
       W_1(l_2)=T12;
    
    end
       end
    end
   
end


for ye3=0.05:0.05:8       %Plane 3
    for ze3=1:0.05:3       %xe3=0
   l_3=l_3+1;
     
      for xe1=0.2:0.2:4       %Plane 1
       for ze1=1:0.2:3          %ye1=0
         l_4=l_4+1;
   
   R71=sqrt((xt-xe3)^2+(yt-ye3)^2+(zt-ze3)^2);
   R73=sqrt((xe3-xe1)^2+(ye3-ye1)^2+(ze3-ze1)^2);
   R74=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R75=R71+R73+R74;
  
  %angel for Second Ref
  theta3_11=acos(abs(zt-ze3)/R71);
  theta3_12=acos(abs(xt-xe3)/R71);              %xt-xe3=2-0=2
  theta3_13=acos(abs(xe1-xe3)/R73);             %% OR ye1-ye3   
  theta3_14=acos(abs(ze3-ze1)/R73);
  theta3_15=acos(abs(zr-ze1)/R74);          
  theta3_16=acos(abs(yr-ye1)/R74);              %yr-ye1=1-0=1
 
  Pres13=(((n+1)*(n+1))^2/(8*(pi)^3*(R71)^2*(R73)^2*(R74)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta3_11))^n*(cos(theta3_12))^n*cos(theta3_13)*(cos(theta3_14))^n*cos(theta3_15)*cos(theta3_16);
   T13=R75/S;
   K_2(l_4)=Pres13;
   W_2(l_4)=T13;
    end
       end
    end  
end

for ye3=0.05:0.05:8       %Plane 3
    for ze3=1:0.05:3       %xe3=0
      l_5=l_5+1;
     
      for xe2=0.2:0.2:4     %Plane 2
        for ze2=1:0.2:3        %ye2=8                                                                                             
           l_6=l_6+1;
   
    R76=sqrt((xt-xe3)^2+(yt-ye3)^2+(zt-ze3)^2);
    R78=sqrt((xe3-xe2)^2+(ye3-ye2)^2+(ze3-ze2)^2);
    R79=sqrt((xr-xe2)^2+(yr-ye2)^2+(zr-ze2)^2);
     R80=R76+R78+R79;
    
      %angel for Second Ref
    theta3_19=acos(abs(xt-xe3)/R76);          %xt-xe3=2-0=2
    theta3_20=acos(abs(zt-ze3)/R76);
    theta3_21=acos(abs(xe2-xe3)/R78);         %% OR ye2-ye3
    theta3_22=acos(abs(ze3-ze2)/R78);
    theta3_23=acos(abs(zr-ze2)/R79); 
    theta3_24=acos(abs(yr-ye2)/R79);          %yr-ye2=8-1=7
    
    Pres14=(((n+1)*(n+1))^2/(8*(pi)^3*(R76)^2*(R78)^2*(R79)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta3_19))^n*(cos(theta3_20))^n*cos(theta3_21)*(cos(theta3_22))^n*cos(theta3_23)*cos(theta3_24);
     T14=R80/S;
     K_3(l_6)=Pres14;
     W_3(l_6)=T14;
    
    end
       end
    end
   
end


for ye3=0.05:0.05:8         %Plane 3
    for ze3=1:0.05:3         %xe3=0
      l_7=l_7+1;
  
       for ye4=0.2:0.2:8        %Plane 4
        for ze4=1:0.2:3           %xe4=4
          l_8=l_8+1;
   
    R81=sqrt((xt-xe3)^2+(yt-ye3)^2+(zt-ze3)^2);
    R83=sqrt((xe3-xe4)^2+(ye3-ye4)^2+(ze3-ze4)^2);
    R84=sqrt((xr-xe4)^2+(yr-ye4)^2+(zr-ze4)^2);
    R85=R81+R83+R84;
  
    %angel for Second Ref
    theta3_27=acos(abs(xt-xe3)/R81);            %xt-xe3=2-0=2
    theta3_28=acos(abs(zt-ze3)/R81);
    theta3_29=acos(abs(ze-ze3)/R83);             
    theta3_30=acos(abs(xe-xe4)/R83);
    theta3_31=acos(abs(zr-ze4)/R84);             
    theta3_32=acos(abs(xr-xe4)/R84);             %xr-xe4=1-4=3 
  
   Pres15=(((n+1)*(n+1))^2/(8*(pi)^3*(R81)^2*(R83)^2*(R84)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta3_27))^n*(cos(theta3_28))^n*cos(theta3_29)*(cos(theta3_30))^n*cos(theta3_31)*cos(theta3_32);
    T15=R85/S; 
    K_4(l_8)=Pres15; 
    W_4(l_8)=T15;
   
     end
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         % ( plane 4 with all )

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ye4=0.05:0.05:8    %Plane 4
    for ze4=1:0.05:3    %xe4=4
   m_1=m_1+1;
   
      for xe=0.2:0.2:4        %ceiling
       for ye=0.2:0.2:8         %ze=3                                                                                               
         m_2=m_2+1;
   
    R86=sqrt((xt-xe4)^2+(yt-ye4)^2+(zt-ze4)^2);
    R88=sqrt((xe-xe4)^2+(ye-ye4)^2+(ze-ze4)^2);
    R89=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R90=R86+R88+R89;
    
    %angel for Second Ref
    theta4_3=acos(abs(xt-xe4)/R86);             %xe4-xt=4-2=2  
    theta4_4=acos(abs(zt-ze4)/R86);
    theta4_5=acos(abs(ze-ze4)/R88);
    theta4_6=acos(abs(xe-xe4)/R88);            %OR ye-ye14
    theta4_7=acos(2/R89);             
    theta4_8=acos(2/R89);     
     
    Pres16=(((n+1)*(n+1))^2/(8*(pi)^3*(R86)^2*(R88)^2*(R89)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta4_3))^n*(cos(theta4_4))^n*cos(theta4_5)*(cos(theta4_6))^n*cos(theta4_7)*cos(theta4_8);
     T16=R90/S;
     Y_1(m_2)=Pres16;
     Z_1(m_2)=T16;
    
    end
       end
    end
   
end

for ye4=0.05:0.05:8    %Plane 4
    for ze4=1:0.05:3    %xe4=4
        m_3=m_3+1;

   for xe1=0.2:0.2:4     %Plane 1
    for ze1=1:0.2:3          %ye1=0
      m_4=m_4+1;
   
   R91=sqrt((xt-xe4)^2+(yt-ye4)^2+(zt-ze4)^2);
   R93=sqrt((xe4-xe1)^2+(ye4-ye1)^2+(ze4-ze1)^2);
   R94=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R95=R91+R93+R94;
  
  %angel for Second Ref
  theta4_11=acos(abs(zt-ze4)/R91);
  theta4_12=acos(abs(xt-xe4)/R91);             %xe4-xt=4-2=2
  theta4_13=acos(abs(xe1-xe4)/R93);            
  theta4_14=acos(abs(ze1-ze4)/R93);
  theta4_15=acos(abs(yr-ye1)/R94);             %yr-ye1=1-0=1
  theta4_16=acos(abs(zr-ze1)/R94);          
  
  Pres17=(((n+1)*(n+1))^2/(8*(pi)^3*(R91)^2*(R93)^2*(R94)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta4_11))^n*(cos(theta4_12))^n*cos(theta4_13)*(cos(theta4_14))^n*cos(theta4_15)*cos(theta4_16);
   T17=R95/S;
   Y_2(m_4)=Pres17;
   Z_2(m_4)=T17;
    end
       end
    end  
end

for ye4=0.05:0.05:8    %Plane 4
    for ze4=1:0.05:3    %xe4=4
   m_5=m_5+1;
  
    for xe2=0.2:0.2:4     %Plane 2
     for ze2=1:0.2:3        %ye2=8                                                                                             
         m_6=m_6+1;
   
   R96=sqrt((xt-xe4)^2+(yt-ye4)^2+(zt-ze4)^2);
   R98=sqrt((xe4-xe2)^2+(ye4-ye2)^2+(ze4-ze2)^2);
   R99=sqrt((xr-xe2)^2+(yr-ye2)^2+(zr-ze2)^2);
    R100=R96+R98+R99;
    
    %angel for Second Ref
   theta4_19=acos(abs(zt-ze2)/R96);
   theta4_20=acos(abs(xt-xe4)/R96);             %xe4-xt=4-2=2
   theta4_21=acos(abs(xe2-xe4)/R98);     
   theta4_22=acos(abs(ze4-ze2)/R98);
   theta4_23=acos(abs(yr-ye2)/R99);             %yr-ye2=8-1=7
   theta4_24=acos(abs(zr-ze2)/R99);   
  
     
   Pres18=(((n+1)*(n+1))^2/(8*(pi)^3*(R96)^2*(R98)^2*(R99)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta4_19))^n*(cos(theta4_20))^n*cos(theta4_21)*(cos(theta4_22))^n*cos(theta4_23)*cos(theta4_24);
    T18=R100/S;
    Y_3(m_6)=Pres18;
    Z_3(m_6)=T18;
    
    end
       end
    end
   
end


   for ye4=0.05:0.05:8    %Plane 4
    for ze4=1:0.05:3    %xe4=4
   m_7=m_7+1;
   
    for ye3=0.2:0.2:8       %Plane 3
     for ze3=1:0.2:3          %xe3=0
       m_8=m_8+1;

     R101=sqrt((xt-xe4)^2+(yt-ye4)^2+(zt-ze4)^2);
     R103=sqrt((xe4-xe3)^2+(ye4-ye3)^2+(ze4-ze3)^2);
     R104=sqrt((xr-xe3)^2+(yr-ye3)^2+(zr-ze3)^2);
     R105=R101+R103+R104;
     
      %angel for Second Ref
      theta4_27=acos(abs(zt-ze4)/R101);
      theta4_28=acos(abs(xt-xe4)/R101);        %xe4-xt=4-2=2
      theta4_29=acos(abs(xe3-xe4)/R103);        
      theta4_30=acos(abs(ze-ze3)/R103);
      theta4_31=acos(abs(zr-ze3)/R104);      
      theta4_32=acos(abs(xr-xe3)/R104);         %xr-xe3=1-0=1 
      
    Pres19=(((n+1)*(n+1))^2/(8*(pi)^3*(R101)^2*(R103)^2*(R104)^2))*Tf*Tc*Pt*dA*dA2*P1*P2*Ar*(cos(theta4_27))^n*(cos(theta4_28))^n*cos(theta4_29)*(cos(theta4_30))^n*cos(theta4_31)*cos(theta4_32);
     T19=R105/S;
     Y_4(m_8)=Pres19; 
     Z_4(m_8)=T19;
 
     end
    end
   end
   end

   
   
   
   B(find(isnan(B)))=0;                  % Ceiling + plane1

  new_plosC1=zeros(1,80);      

    for q=1:80  

      cc=(floor((C*1e9)))-q;

       [rows,cols]=find(cc==0);

       if(cols>0)

       prq=sum(B(cols));

       new_plosC1(q)=new_plosC1(q)+prq;

       end            

    end
    
      
   B_1(find(isnan(B_1)))=0;                  % Ceiling + plane2

  new_plosC2=zeros(1,80);      

    for q=1:80  

      cc1=(floor((C_1*1e9)))-q;

       [rows,cols]=find(cc1==0);

       if(cols>0)

       prq=sum(B_1(cols));

       new_plosC2(q)=new_plosC2(q)+prq;

       end            

    end

     B_2(find(isnan(B_2)))=0;                  % Ceiling + plane3

  new_plosC3=zeros(1,80);      

    for q=1:80  

      cc2=(floor((C_2*1e9)))-q;

       [rows,cols]=find(cc2==0);

       if(cols>0)

       prq=sum(B_2(cols));

       new_plosC3(q)=new_plosC3(q)+prq;

       end            

    end
    
       B_3(find(isnan(B_3)))=0;                  % Ceiling + plane4

  new_plosC4=zeros(1,80);      

    for q=1:80  

      cc3=(floor((C_3*1e9)))-q;

       [rows,cols]=find(cc3==0);

       if(cols>0)

       prq=sum(B_3(cols));

       new_plosC4(q)=new_plosC4(q)+prq;

       end            

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     F_1(find(isnan(F_1)))=0;                  % plane1 + Ceiling

   new_plosP1=zeros(1,80);      

    for q=1:80  

      gg1=(floor((G_1*1e9)))-q;

       [rows,cols]=find(gg1==0);

       if(cols>0)

       prq=sum(F_1(cols));

       new_plosP1(q)=new_plosP1(q)+prq;

       end            

    end
    
     F_2(find(isnan(F_2)))=0;                  % plane1 + plane2

   new_plosP2=zeros(1,80);      

    for q=1:80  

      gg2=(floor((G_2*1e9)))-q;

       [rows,cols]=find(gg2==0);

       if(cols>0)

       prq=sum(F_2(cols));

       new_plosP2(q)=new_plosP2(q)+prq;

       end            

    end

    F_3(find(isnan(F_3)))=0;                  % plane1 + plane3

   new_plosP3=zeros(1,80);      

    for q=1:80  

      gg3=(floor((G_3*1e9)))-q;

       [rows,cols]=find(gg3==0);

       if(cols>0)

       prq=sum(F_3(cols));

       new_plosP3(q)=new_plosP3(q)+prq;

       end            

    end
    
    F_4(find(isnan(F_4)))=0;                  % plane1 + plane4

   new_plosP4=zeros(1,80);      

    for q=1:80  

      gg4=(floor((G_4*1e9)))-q;

       [rows,cols]=find(gg4==0);

       if(cols>0)

       prq=sum(F_4(cols));

       new_plosP4(q)=new_plosP4(q)+prq;

       end            

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     I_1(find(isnan(I_1)))=0;                  % plane2 + Ceiling

   new_plosP5=zeros(1,80);      

    for q=1:80  

      hh1=(floor((H_1*1e9)))-q;

       [rows,cols]=find(hh1==0);

       if(cols>0)

       prq=sum(I_1(cols));

       new_plosP5(q)=new_plosP5(q)+prq;

       end            

    end
    
    I_2(find(isnan(I_2)))=0;                  % plane2 + plane1

   new_plosP6=zeros(1,80);      

    for q=1:80  

      hh2=(floor((H_2*1e9)))-q;

       [rows,cols]=find(hh2==0);

       if(cols>0)

       prq=sum(I_2(cols));

       new_plosP6(q)=new_plosP6(q)+prq;

       end            

    end

    I_3(find(isnan(I_3)))=0;                  % plane2 + plane3

   new_plosP7=zeros(1,80);      

    for q=1:80  

      hh3=(floor((H_3*1e9)))-q;

       [rows,cols]=find(hh3==0);

       if(cols>0)

       prq=sum(I_3(cols));

       new_plosP7(q)=new_plosP7(q)+prq;

       end            

    end

     I_4(find(isnan(I_4)))=0;                  % plane2 + plane4

   new_plosP8=zeros(1,80);      

    for q=1:80  

      hh4=(floor((H_4*1e9)))-q;

       [rows,cols]=find(hh4==0);

       if(cols>0)

       prq=sum(I_4(cols));

       new_plosP8(q)=new_plosP8(q)+prq;

       end            

    end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     K_1(find(isnan(K_1)))=0;                  % plane3 + Ceiling

   new_plosP9=zeros(1,80);      

    for q=1:80  

      ww1=(floor((W_1*1e9)))-q;

       [rows,cols]=find(ww1==0);

       if(cols>0)

       prq=sum(K_1(cols));

       new_plosP9(q)=new_plosP9(q)+prq;

       end            

    end
    
       K_2(find(isnan(K_2)))=0;                  % plane3 + plane1

   new_plosP10=zeros(1,80);      

    for q=1:80  

      ww2=(floor((W_2*1e9)))-q;

       [rows,cols]=find(ww2==0);

       if(cols>0)

       prq=sum(K_2(cols));

       new_plosP10(q)=new_plosP10(q)+prq;

       end            

    end
    
       K_3(find(isnan(K_3)))=0;                  % plane3 + plane2

   new_plosP11=zeros(1,80);      

    for q=1:80  

      ww3=(floor((W_3*1e9)))-q;

       [rows,cols]=find(ww3==0);

       if(cols>0)

       prq=sum(K_3(cols));

       new_plosP11(q)=new_plosP11(q)+prq;

       end            

    end
    
       K_4(find(isnan(K_4)))=0;                  % plane3 + plane4

   new_plosP12=zeros(1,80);      

    for q=1:80  

      ww4=(floor((W_4*1e9)))-q;

       [rows,cols]=find(ww4==0);

       if(cols>0)

       prq=sum(K_4(cols));

       new_plosP12(q)=new_plosP12(q)+prq;

       end            

    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     Y_1(find(isnan(Y_1)))=0;                  % plane4 + Ceiling

   new_plosP13=zeros(1,80);      

    for q=1:80  

      zz1=(floor((Z_1*1e9)))-q;

       [rows,cols]=find(zz1==0);

       if(cols>0)

       prq=sum(Y_1(cols));

       new_plosP13(q)=new_plosP13(q)+prq;

       end            

    end
    
    Y_2(find(isnan(Y_2)))=0;                  % plane4 + plane1

   new_plosP14=zeros(1,80);      

    for q=1:80  

      zz2=(floor((Z_2*1e9)))-q;

       [rows,cols]=find(zz2==0);

       if(cols>0)

       prq=sum(Y_2(cols));

       new_plosP14(q)=new_plosP14(q)+prq;

       end            

    end
    
      Y_3(find(isnan(Y_3)))=0;                  % plane4 + plane2

   new_plosP15=zeros(1,80);      

    for q=1:80  

      zz3=(floor((Z_3*1e9)))-q;

       [rows,cols]=find(zz3==0);

       if(cols>0)

       prq=sum(Y_3(cols));

       new_plosP15(q)=new_plosP15(q)+prq;

       end            

    end
    
     Y_4(find(isnan(Y_4)))=0;                  % plane4 + plane3

   new_plosP16=zeros(1,80);      

    for q=1:80  

      zz4=(floor((Z_4*1e9)))-q;

       [rows,cols]=find(zz4==0);

       if(cols>0)

       prq=sum(Y_4(cols));

       new_plosP16(q)=new_plosP16(q)+prq;

       end            

    end
    
SecondReflection_plos=new_plosC1+new_plosC2+new_plosC3+new_plosC4+new_plosP1+new_plosP2+new_plosP3+new_plosP4+new_plosP5+new_plosP6+new_plosP7+new_plosP8+new_plosP9+new_plosP10+new_plosP11+new_plosP12+new_plosP13+new_plosP14+new_plosP15+new_plosP16;
     
     Total_First_Second=FirstReflection_plos+SecondReflection_plos;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
                  %%%%%%%%%% Delay Spread  %%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     Optical_Power=FirstReflection_plos+SecondReflection_plos;
    
       Time=1:80;
mean_delay=(sum((Time.*(Optical_Power).^2)))/(sum((Optical_Power).^2));
Delay=sqrt(sum((((Time-mean_delay).^2).*(Optical_Power).^2))/sum((Optical_Power).^2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
