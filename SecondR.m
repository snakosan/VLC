
clear all;
FOV=90*pi/180; % field of view of the RX
xr=1; yr=1; zr=1;
xt=2; yt=4; zt=1;
n=1;  % radiation beam
Tf=1; % Filter Ttansmission
Tc=1; % Concentrator Ttansmission
Pt=1e-3; % optical power
Ar=1e-4; % collection area
dA= (5*5)*10^(-4);
dA2= (20*20)*10^(-4);
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
A_3=zeros(1,12800);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ( Ceiling with all Planes )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 

for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3                                                                                               
   i=i+1;
   
  R1=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R2=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  

  %angel for First Ref
  theta1=acos(2/R1);
  phi1=theta1;
  theta2=acos(2/R2);
  phi2=theta2;
  
  
  Pre=((n+1)*(n+1)/(4*(pi)^2*(R1)^2*(R2)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta1))^n*cos(phi1)*(cos(theta2))^n*cos(phi2);
  
   for xe1=0.05:0.05:4     %Plane 1
    for ze1=1:0.05:3        %ye1=0
   i_1=i_1+1;
   
   R1=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R3=sqrt((xe-xe1)^2+(ye-ye1)^2+(ze-ze1)^2);
   R4=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R=R1+R2+R3+R4;
  
  %angel for Second Ref
  theta11=acos(2/R1);
  theta12=acos(2/R1);            
  theta13=acos(2/R3);             
  theta14=acos(abs(ze-ze1)/R3);   
  theta15=acos(2/R4);             %zr-ze=1-3=2
  theta16=acos(abs(ze-ze1)/R4); 
  
  Pres=(((n+1)*(n+1))^2/(8*(pi)^3*(R1)^2*(R3)^2*(R4)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta11))^n*(cos(theta12))^n*cos(theta13)*(cos(theta14))^n*cos(theta15)*cos(theta16);
   T=R/S;
   A(i)=Pre;
   B(i_1)=Pres;
   C(i_1)=T;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3
   i_2=i_2+1;
   
 
  R11=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R12=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  

  %angel for First Ref
  theta3=acos(2/R11);
  phi3=theta3;
  theta4=acos(2/R12);
  phi4=theta4;
  

  Pre1=((n+1)*(n+1)/(4*(pi)^2*(R11)^2*(R12)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3))^n*cos(phi3)*(cos(theta4))^n*cos(phi4);
  
  
   for xe12=0.05:0.05:4     %Plane 2
    for ze12=1:0.05:3        %ye2=8
   i_3=i_3+1;
   
   R11=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R13=sqrt((xe-xe12)^2+(ye-ye2)^2+(ze-ze12)^2);
   R14=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
   R15=R11+R12+R13+R14;
  
   %angel for Second Ref
  theta17=acos(2/R11);
  theta18=acos(2/R11);             
  theta19=acos(2/R13);             
  theta20=acos(abs(ze-ze12)/R13);   
  theta21=acos(2/R14);              
  theta22=acos(abs(ze-ze12)/R14);              
  
  Pres1=(((n+1)*(n+1))^2/(8*(pi)^3*(R11)^2*(R13)^2*(R14)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta17))^n*(cos(theta18))^n*cos(theta19)*(cos(theta20))^n*cos(theta21)*cos(theta22);
    T1=R15/S;
   A_1(i_2)=Pre1;
   B_1(i_3)=Pres1; 
   C_1(i_3)=T1;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3
   i_4=i_4+1;
   
  R16=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R17=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  

  %angel for First Ref
  theta5=acos(2/R16);
  phi5=theta5;
  theta6=acos(2/R17);
  phi6=theta6;
  
 
  Pre2=((n+1)*(n+1)/(4*(pi)^2*(R16)^2*(R17)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta5))^n*cos(phi5)*(cos(theta6))^n*cos(phi6);
  
  
   for ye13=0.05:0.05:8    %Plane 3
     for ze13=1:0.05:3    %xe3=0
           i_5=i_5+1;

     R16=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
     R18=sqrt((xe-xe3)^2+(ye-ye13)^2+(ze-ze13)^2);
     R19=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
     R20=R16+R17+R18+R19;

        %angel for Second Ref
      theta23=acos(2/R16);
      theta24=acos(2/R16);             %xe3-xt=2-0=2
      theta25=acos(2/R18);             %xe3-xt=2-0=2
      theta26=acos(abs(ze-ze13)/R18);
      theta27=acos(2/R19);             %zr-ze=1-3=2
      theta28=acos(abs(ze-ze13)/R19);           

    Pres2=(((n+1)*(n+1))^2/(8*(pi)^3*(R16)^2*(R18)^2*(R19)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta23))^n*(cos(theta24))^n*cos(theta25)*(cos(theta26))^n*cos(theta27)*cos(theta28);
         T2=R20/S;
         A_2(i_4)=Pre2;
         B_2(i_5)=Pres2; 
         C_2(i_5)=T2;
    end
       end
    end
   
end


for xe=0.05:0.05:4      %ceiling
    for ye=0.05:0.05:8    %ze=3                                                                                               
   i_6=i_6+1;
   
  R21=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
  R22=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
  

  %angel for First Ref
  theta7=acos(2/R21);
  phi7=theta7;
  theta8=acos(2/R22);
  phi8=theta8;
  
  
  Pre3=((n+1)*(n+1)/(4*(pi)^2*(R21)^2*(R22)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta7))^n*cos(phi7)*(cos(theta8))^n*cos(phi8);
  A_3(i_6)=Pre3;
  
   
   for ye14=0.05:0.05:8    %Plane 4
     for ze14=1:0.05:3    %xe4=4
   i_7=i_7+1;
   
   R21=sqrt((xt-xe)^2+(yt-ye)^2+(zt-ze)^2);
   R23=sqrt((xe-xe4)^2+(ye-ye14)^2+(ze-ze14)^2);
   R24=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);
   R25=R21+R22+R23+R24;
  
    %angel for Second Ref
  theta29=acos(2/R21);
  theta30=acos(2/R21);             %xe4-xt=4-2=2
  theta31=acos(2/R23);             %xe4-xt=4-2=2
  theta32=acos(abs(ze-ze14)/R23);
  theta33=acos(2/R24);             %zr-ze=1-3=2
  theta34=acos(abs(ze-ze14)/R24);             
  
  Pres3=(((n+1)*(n+1))^2/(8*(pi)^3*(R21)^2*(R23)^2*(R24)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta29))^n*(cos(theta30))^n*cos(theta31)*(cos(theta32))^n*cos(theta33)*cos(theta34);
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
   
    R26=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R27=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
    
    
    %angel for First Ref
    theta35=acos(4/R26);            %yt-ye1=4-0=4
    phi9=acos(abs(zt-ze1)/R26);
    theta36=acos(1/R27);            %yr-ye1=1-0=1
    phi10=acos(abs(zr-ze1)/R27);
    
    Pre4=((n+1)*(n+1)/(4*(pi)^2*(R26)^2*(R27)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta35))^n*cos(phi9)*(cos(theta36))^n*cos(phi10);
    E_1(j_1)=Pre4;
    
    
   for xe=0.05:0.05:4      %ceiling
   for ye=0.05:0.05:8      %ze=3                                                                                               
   j_2=j_2+1;
   
   R26=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
   R28=sqrt((xe-xe1)^2+(ye-ye1)^2+(ze-ze1)^2);
   R29=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R30=R26+R27+R28+R29;
    
     %angel for Second Ref
  theta37=acos(abs(zt-ze1)/R26);
  theta38=acos(4/R26);             %yt-ye1=4-0=4
  theta39=acos(abs(ze-ze1)/R28);
  theta40=acos(abs(ye-ye1)/R28);    %ye1-ye
  theta41=acos(2/R29);            
  theta42=acos(2/R29);     
     
  Pres4=(((n+1)*(n+1))^2/(8*(pi)^3*(R26)^2*(R28)^2*(R29)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta37))^n*(cos(theta38))^n*cos(theta39)*(cos(theta40))^n*cos(theta41)*cos(theta42);
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
   
    R31=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R32=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
    
    
    %angel for First Ref
    theta43=acos(4/R31);
    phi11=acos(abs(zt-ze1)/R31);
    theta44=acos(1/R32);
    phi12=acos(abs(zr-ze1)/R32);
    
    Pre5=((n+1)*(n+1)/(4*(pi)^2*(R31)^2*(R32)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta43))^n*cos(phi11)*(cos(theta44))^n*cos(phi12);
    E_2(j_3)=Pre5;
    
    
   for xe12=0.05:0.05:4     %Plane 2
    for ze12=1:0.05:3        %ye2=8                                                                                             
   j_4=j_4+1;
   
   R31=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
   R33=sqrt((xe1-xe12)^2+(ye1-ye2)^2+(ze1-ze12)^2);
   R34=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
   R35=R31+R32+R33+R34;
    
   %angel for Second Ref
  theta45=acos(abs(zt-ze1)/R31);
  theta46=acos(4/R31);                   %ye1-yt=4-0=4
  theta47=acos(8/R33);                  %ye2-ye1=8-0=8
  theta48=acos(abs(ze1-ze12)/R33);
  theta49=acos(abs(zr-ze12)/R34);           
  theta50=acos(7/R34);                   %ye2-yr=8-1=7    
  
     
  Pres5=(((n+1)*(n+1))^2/(8*(pi)^3*(R31)^2*(R33)^2*(R34)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta45))^n*(cos(theta46))^n*cos(theta47)*(cos(theta48))^n*cos(theta49)*cos(theta50);
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
   
    R36=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R37=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
    
    
    %angel for First Ref
    theta51=acos(4/R36);          %yt-ye1=4-0=4
    phi13=acos(abs(zt-ze1)/R36);
    theta52=acos(1/R37);          %yr-ye1=1-0=1
    phi14=acos(abs(zr-ze1)/R37);
    
    Pre6=((n+1)*(n+1)/(4*(pi)^2*(R36)^2*(R37)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta51))^n*cos(phi13)*(cos(theta52))^n*cos(phi14);
    E_3(j_5)=Pre6;
    
    
     for ye13=0.05:0.05:8    %Plane 3
     for ze13=1:0.05:3    %xe3=0
           j_6=j_6+1;

      R36=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
      R38=sqrt((xe1-xe3)^2+(ye1-ye13)^2+(ze1-ze13)^2);
      R39=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
      R40=R36+R37+R38+R39;

      %angel for Second Ref
      theta53=acos(abs(zt-ze1)/R36);
      theta54=acos(4/R36);             %yt-ye1=4-0=4
      theta55=acos(abs(xe1-xe3)/R38);        
      theta56=acos(abs(ze1-ze13)/R38);
      theta57=acos(abs(zr-ze13)/R39);            
      theta58=acos(abs(yr-ye13)/R39);         

     Pres6=(((n+1)*(n+1))^2/(8*(pi)^3*(R36)^2*(R38)^2*(R39)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta53))^n*(cos(theta54))^n*cos(theta55)*(cos(theta56))^n*cos(theta57)*cos(theta58);
        T6=R40/S;
        E_3(j_5)=Pre6;
        F_3(j_6)=Pres6; 
        G_3(j_6)=T6;
    end
       end
    end
   
end

   for xe1=0.05:0.05:4   %Plane 1
    for ze1=1:0.05:3      %ye1=0
   j_7=j_7+1;
   
    R41=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R42=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
    
    
    %angel for First Ref
    theta59=acos(4/R41);         %yt-ye1=4-0=4
    phi15=acos(abs(zt-ze1)/R41);
    theta60=acos(1/R42);         %yr-ye1=1-0=1
    phi16=acos(abs(zr-ze1)/R42);
    
    Pre7=((n+1)*(n+1)/(4*(pi)^2*(R41)^2*(R42)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta59))^n*cos(phi15)*(cos(theta60))^n*cos(phi16);
    E_4(j_7)=Pre7;
    
     
       for ye14=0.05:0.05:8    %Plane 4
          for ze14=1:0.05:3    %xe4=4
             j_8=j_8+1;
   
    R41=sqrt((xt-xe1)^2+(yt-ye1)^2+(zt-ze1)^2);
    R43=sqrt((xe1-xe4)^2+(ye1-ye14)^2+(ze1-ze14)^2);
    R44=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);
    R45=R41+R42+R43+R44;
   
       %angel for Second Ref
    theta61=acos(abs(zt-ze1)/R41);
    theta62=acos(4/R41);             %yt-ye1=4-0=4
    theta63=acos(abs(xe1-xe4)/R43);         
    theta64=acos(abs(ze1-ze14)/R43);
    theta65=acos(abs(zr-ze14)/R44);             
    theta66=acos(abs(yr-ye14)/R44);             
  
  Pres7=(((n+1)*(n+1))^2/(8*(pi)^3*(R41)^2*(R43)^2*(R44)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta61))^n*(cos(theta62))^n*cos(theta63)*(cos(theta64))^n*cos(theta65)*cos(theta66);
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
    
    
    for xe12=0.05:0.05:4        %Plane 2
    for ze12=1:0.05:3        %ye2=8
      k_1=k_1+1;
   
    R46=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
    R47=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
   
    %angel for First Ref
    theta67=acos(6/46);      %ye2-yt=8-4=6    
    phi17=acos(abs(zt-ze12)/R46);
    theta68=acos(7/R47);     %ye2-yr=8-1=7 
    phi18=acos(abs(zr-ze12)/R47);
    
    Pre8=((n+1)*(n+1)/(4*(pi)^2*(R46)^2*(R47)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta67))^n*cos(phi17)*(cos(theta68))^n*cos(phi18);
     D_1(k_1)=Pre8;
    
    
    for xe=0.05:0.05:4       %ceiling
   for ye=0.05:0.05:8        %ze=3                                                                                               
     k_2=k_2+1;
     
   R46=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
   R48=sqrt((xe12-xe)^2+(ye2-ye)^2+(ze12-ze)^2);
   R49=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R50=R46+R47+R48+R49;
    
     %angel for Second Ref
  theta69=acos(6/46);      %ye2-yt=8-4=6
  theta70=acos(abs(zt-ze12)/R46);
  theta71=acos(abs(ze-ze12)/R48);
  theta72=acos(abs(ye-ye2)/R48);       %ye2-ye 
  theta73=acos(2/R49);       
  theta74=acos(2/R49);     
     
  Pres8=(((n+1)*(n+1))^2/(8*(pi)^3*(R46)^2*(R48)^2*(R49)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta69))^n*(cos(theta70))^n*cos(theta71)*(cos(theta72))^n*cos(theta73)*cos(theta74);
    T8=R50/S;
    I_1(k_2)=Pres8;
    H_1(k_2)=T8;
    
    end
       end
    end
   
end
    
for xe12=0.05:0.05:4          %Plane 2
    for ze12=1:0.05:3           %ye2=8
      k_3=k_3+1;
   
    R51=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
    R52=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
   
    %angel for First Ref
    theta75=acos(4/51);        %ye2-yt=8-4=4
    phi19=acos(abs(zt-ze12)/R51);
    theta76=acos(7/R52);       %ye2-yr=8-1=7
    phi20=acos(abs(zr-ze12)/R52);
    
    Pre9=((n+1)*(n+1)/(4*(pi)^2*(R51)^2*(R52)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta75))^n*cos(phi19)*(cos(theta76))^n*cos(phi20);
     D_2(k_3)=Pre9;
      
      for xe1=0.05:0.05:4       %Plane 1
    for ze1=1:0.05:3              %ye1=0
      k_4=k_4+1;
   
   R51=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
   R53=sqrt((xe12-xe1)^2+(ye2-ye1)^2+(ze12-ze1)^2);
   R54=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R55=R51+R52+R53+R54;
  
  %angel for Second Ref
  theta77=acos(4/51);            %ye2-yt=8-4=4
  theta78=acos(abs(zt-ze12)/R51);
  theta79=acos(8/R53);             %ye1-y2=8-0=8
  theta80=acos(abs(ze12-ze1)/R53);
  theta81=acos(1/R54);             %yr-ye1=1-0=1
  theta82=acos(abs(zr-ze1)/R54);       
  
  Pres9=(((n+1)*(n+1))^2/(8*(pi)^3*(R51)^2*(R53)^2*(R54)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta77))^n*(cos(theta78))^n*cos(theta79)*(cos(theta80))^n*cos(theta81)*cos(theta82);
   T9=R55/S;
   I_2(k_4)=Pres9;
   H_2(k_4)=T9;
    end
       end
    end  
end

for xe12=0.05:0.05:4          %Plane 2
 for ze12=1:0.05:3             %ye2=8
  k_5=k_5+1;

   R56=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
   R57=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
    
   %angel for First Ref
    theta83=acos(4/56);        %ye2-yt=8-4=4
    phi21=acos(abs(zt-ze12)/R56);
    theta84=acos(7/R57);       %ye2-yr=8-1=7
    phi22=acos(abs(zr-ze12)/R57);

   Pre10=((n+1)*(n+1)/(4*(pi)^2*(R56)^2*(R57)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta83))^n*cos(phi21)*(cos(theta84))^n*cos(phi22);
     D_3(k_5)=Pre10;

     
    for ye13=0.05:0.05:8       %Plane 3
     for ze13=1:0.05:3          %xe3=0
      k_6=j_6+1;

     R56=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
     R58=sqrt((xe12-xe3)^2+(ye2-ye13)^2+(ze12-ze13)^2);
     R59=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
     R60=R56+R57+R58+R59;
     
       %angel for Second Ref
      theta85=acos(4/56);           %ye2-yt=8-4=4
      theta86=acos(abs(zt-ze12)/R56);
      theta87=acos(abs(ye2-ye13)/R58);      
      theta88=acos(abs(ze12-ze13)/R58);
      theta89=acos(abs(zr-ze13)/R59);   
      theta90=acos(abs(yr-ye13)/R59);            
      
      Pres10=(((n+1)*(n+1))^2/(8*(pi)^3*(R56)^2*(R58)^2*(R59)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta85))^n*(cos(theta86))^n*cos(theta87)*(cos(theta88))^n*cos(theta89)*cos(theta90);
  
          T10=R60/S;
          I_3(k_6)=Pres10; 
          H_3(k_6)=T10;
 
     end
    end
   end
end
 

for xe12=0.05:0.05:4         %Plane 2
 for ze12=1:0.05:3            %ye2=8
  k_7=k_7+1;

   R61=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
   R62=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
    
   %angel for First Ref
    theta91=acos(4/61);        %ye2-yt=8-4=4
    phi23=acos(abs(zt-ze12)/R61);
    theta92=acos(7/R62);       %ye2-yr=8-1=7
    phi24=acos(abs(zr-ze12)/R62);
    
    Pre11=((n+1)*(n+1)/(4*(pi)^2*(R61)^2*(R62)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta91))^n*cos(phi23)*(cos(theta92))^n*cos(phi24);
     D_4(k_7)=Pre11;
     
       for ye14=0.05:0.05:8        %Plane 4
        for ze14=1:0.05:3            %xe4=4
          k_8=k_8+1;
   
    R61=sqrt((xt-xe12)^2+(yt-ye2)^2+(zt-ze12)^2);
    R63=sqrt((xe12-xe4)^2+(ye2-ye14)^2+(ze12-ze14)^2);
    R64=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);
    R65=R61+R62+R63+R64;
  
      %angel for Second Ref
    theta93=acos(4/61);         %ye2-yt=8-4=4
    theta94=acos(abs(zt-ze12)/R61);
    theta95=acos(abs(xe12-xe4)/R63);             
    theta96=acos(abs(ze12-ze14)/R63);
    theta97=acos(abs(zr-ze14)/R64);    
    theta98=acos(abs(yr-ye14)/R64);              
  
    Pres11=(((n+1)*(n+1))^2/(8*(pi)^3*(R61)^2*(R63)^2*(R64)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta93))^n*(cos(theta94))^n*cos(theta95)*(cos(theta96))^n*cos(theta97)*cos(theta98);
     
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

 
 
for ye13=0.05:0.05:8       %Plane 3
    for ze13=1:0.05:3       %xe3=0
   l_1=l_1+1;
   
    R66=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R67=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
  
    %angel for First Ref
    theta3_1=acos(2/R66);      %xt-xe3=2-0=2  
    phi25=acos(abs(zt-ze13)/R66);
    theta3_2=acos(1/R67);      %xr-xe3=1-0=1 
    phi26=acos(abs(zr-ze13)/R67);
    
    Pre12=((n+1)*(n+1)/(4*(pi)^2*(R66)^2*(R67)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3_1))^n*cos(phi25)*(cos(theta3_2))^n*cos(phi26);
     N_1(l_1)=Pre12;
 
    
     for xe=0.05:0.05:4        %ceiling
      for ye=0.05:0.05:8         %ze=3                                                                                               
        l_2=l_2+1;
   
     R66=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
     R68=sqrt((xe-xe3)^2+(ye-ye13)^2+(ze-ze13)^2);
     R69=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
        R70=R66+R67+R68+R69;
    
        %angel for Second Ref
     theta3_3=acos(2/R66);             %xt-xe3=2-0=2  
     theta3_4=acos(abs(zt-ze13)/R66);
     theta3_5=acos(abs(ze-ze13)/R68);
     theta3_6=acos(xe/R68);             % OR ye-ye13  
     theta3_7=acos(2/R69);             
     theta3_8=acos(2/R69);     
     
      Pres12=(((n+1)*(n+1))^2/(8*(pi)^3*(R66)^2*(R68)^2*(R69)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta3_3))^n*(cos(theta3_4))^n*cos(theta3_5)*(cos(theta3_6))^n*cos(theta3_7)*cos(theta3_8);
   
       T12=R70/S;
       K_1(l_2)=Pres12;
       W_1(l_2)=T12;
    
    end
       end
    end
   
end


for ye13=0.05:0.05:8       %Plane 3
    for ze13=1:0.05:3       %xe3=0
   l_3=l_3+1;
   
    R71=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R72=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
  
    %angel for First Ref
    theta3_9=acos(2/R71);            %xt-xe3=2-0=2
    phi27=acos(abs(zt-ze13)/R71);
    theta3_10=acos(1/R72);           %xr-xe3=1-0=1
    phi28=acos(abs(zr-ze13)/R72);
    
    Pre13=((n+1)*(n+1)/(4*(pi)^2*(R71)^2*(R72)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3_9))^n*cos(phi27)*(cos(theta3_10))^n*cos(phi28);
     N_2(l_3)=Pre13;
     
     
      for xe1=0.05:0.05:4       %Plane 1
       for ze1=1:0.05:3          %ye1=0
         l_4=l_4+1;
   
   R71=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
   R73=sqrt((xe3-xe1)^2+(ye13-ye1)^2+(ze13-ze1)^2);
   R74=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R75=R71+R72+R73+R74;
  
  %angel for Second Ref
  theta3_11=acos(abs(zt-ze13)/R71);
  theta3_12=acos(2/R71);                     %xt-xe3=2-0=2
  theta3_13=acos(abs(xe1-xe3)/R73);                    
  theta3_14=acos(abs(ze13-ze1)/R73);
  theta3_15=acos(abs(zr-ze1)/R74);          
  theta3_16=acos(1/R74);                      %yr-ye1=1-0=1
  
  Pres13=(((n+1)*(n+1))^2/(8*(pi)^3*(R71)^2*(R73)^2*(R74)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta3_11))^n*(cos(theta3_12))^n*cos(theta3_13)*(cos(theta3_14))^n*cos(theta3_15)*cos(theta3_16);
   T13=R75/S;
   K_2(l_4)=Pres13;
   W_2(l_4)=T13;
    end
       end
    end  
end

for ye13=0.05:0.05:8       %Plane 3
    for ze13=1:0.05:3       %xe3=0
      l_5=l_5+1;
   
    R76=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R77=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
  
    %angel for First Ref
    theta3_17=acos(2/R76);            %xt-xe3=2-0=2
    phi29=acos(abs(zt-ze13)/R76);
    theta3_18=acos(1/R77);            %xr-xe3=1-0=1
    phi30=acos(abs(zr-ze13)/R77);
    
    Pre14=((n+1)*(n+1)/(4*(pi)^2*(R76)^2*(R77)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3_17))^n*cos(phi29)*(cos(theta3_18))^n*cos(phi30);
     N_3(l_5)=Pre14;
     
     
      for xe12=0.05:0.05:4     %Plane 2
        for ze12=1:0.05:3        %ye2=8                                                                                             
           l_6=l_6+1;
   
    R76=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R78=sqrt((xe3-xe12)^2+(ye13-ye2)^2+(ze13-ze12)^2);
    R79=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
     R80=R76+R77+R78+R79;
    
      %angel for Second Ref
    theta3_19=acos(2/R76);             %xt-xe3=2-0=2
    theta3_20=acos(abs(zt-ze13)/R76);
    theta3_21=acos(abs(xe12-xe3)/R78);
    theta3_22=acos(abs(ze13-ze12)/R78);
    theta3_23=acos(abs(zr-ze12)/R79); 
    theta3_24=acos(7/R79);             %yr-ye2=8-1=7
    
    Pres14=(((n+1)*(n+1))^2/(8*(pi)^3*(R76)^2*(R78)^2*(R79)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta3_19))^n*(cos(theta3_20))^n*cos(theta3_21)*(cos(theta3_22))^n*cos(theta3_23)*cos(theta3_24);
     T14=R80/S;
     K_3(l_6)=Pres14;
     W_3(l_6)=T14;
    
    end
       end
    end
   
end


for ye13=0.05:0.05:8         %Plane 3
    for ze13=1:0.05:3         %xe3=0
      l_7=l_7+1;
   
    R81=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R82=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
  
    %angel for First Ref
    theta3_25=acos(2/R81);            %xt-xe3=2-0=2
    phi31=acos(abs(zt-ze13)/R81);
    theta3_26=acos(1/R82);            %xr-xe3=1-0=1
    phi32=acos(abs(zr-ze13)/R82);
    
    Pre15=((n+1)*(n+1)/(4*(pi)^2*(R81)^2*(R82)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta3_25))^n*cos(phi31)*(cos(theta3_26))^n*cos(phi32);
     N_4(l_7)=Pre15;
     
       for ye14=0.05:0.05:8        %Plane 4
        for ze14=1:0.05:3           %xe4=4
          l_8=l_8+1;
   
    R81=sqrt((xt-xe3)^2+(yt-ye13)^2+(zt-ze13)^2);
    R83=sqrt((xe3-xe4)^2+(ye13-ye14)^2+(ze13-ze14)^2);
    R84=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);
    R85=R81+R82+R83+R84;
  
    %angel for Second Ref
    theta3_27=acos(2/R81);            %xt-xe3=2-0=2
    theta3_28=acos(abs(zt-ze13)/R81);
    theta3_29=acos(abs(ze-ze14)/R83);             
    theta3_30=acos(abs(ze-ze13)/R83);
    theta3_31=acos(abs(zr-ze14)/R84);             
    theta3_32=acos(3/R84);             %xr-xe4=1-4=3 
  
   Pres15=(((n+1)*(n+1))^2/(8*(pi)^3*(R81)^2*(R83)^2*(R84)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta3_27))^n*(cos(theta3_28))^n*cos(theta3_29)*(cos(theta3_30))^n*cos(theta3_31)*cos(theta3_32);
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



for ye14=0.05:0.05:8    %Plane 4
    for ze14=1:0.05:3    %xe4=4
   m_1=m_1+1;
   
    R86=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
    R87=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);

    %angel for First Ref
    theta4_1=acos(2/R86);          %xe4-xt=4-2=2  
    phi33=acos(abs(zt-ze14)/R86);
    theta4_2=acos(3/R87);          %xe4-xr=4-1=3
    phi34=acos(abs(zr-ze14)/R87);
    
    Pre16=((n+1)*(n+1)/(4*(pi)^2*(R86)^2*(R87)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta4_1))^n*cos(phi33)*(cos(theta4_2))^n*cos(phi34);
    X_1(m_1)=Pre16;
  

      for xe=0.05:0.05:4        %ceiling
       for ye=0.05:0.05:8         %ze=3                                                                                               
         m_2=m_2+1;
   
    R86=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
    R88=sqrt((xe-xe4)^2+(ye-ye14)^2+(ze-ze14)^2);
    R89=sqrt((xr-xe)^2+(yr-ye)^2+(zr-ze)^2);
    R90=R86+R87+R88+R89;
    
    %angel for Second Ref
    theta4_3=acos(2/R86);              %xe4-xt=4-2=2  
    theta4_4=acos(abs(zt-ze14)/R86);
    theta4_5=acos(abs(ze-ze14)/R88);
    theta4_6=acos(abs(xe-xe4)/R88);    %OR ye-ye14
    theta4_7=acos(2/R89);             
    theta4_8=acos(2/R89);     
     
    Pres16=(((n+1)*(n+1))^2/(8*(pi)^3*(R86)^2*(R88)^2*(R89)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta4_3))^n*(cos(theta4_4))^n*cos(theta4_5)*(cos(theta4_6))^n*cos(theta4_7)*cos(theta4_8);
     T16=R90/S;
     Y_1(m_2)=Pres16;
     Z_1(m_2)=T16;
    
    end
       end
    end
   
end

for ye14=0.05:0.05:8    %Plane 4
    for ze14=1:0.05:3    %xe4=4
        m_3=m_3+1;
   
    R91=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
    R92=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);

    %angel for First Ref
    theta4_9=acos(2/R91);           %xe4-xt=4-2=2
    phi35=acos(abs(zt-ze14)/R91);
    theta4_10=acos(3/R92);          %xe4-xr=4-1=3
    phi36=acos(abs(zr-ze14)/R92);
    
    Pre17=((n+1)*(n+1)/(4*(pi)^2*(R91)^2*(R92)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta4_9))^n*cos(phi35)*(cos(theta4_10))^n*cos(phi36);
    X_2(m_3)=Pre17;
  

   for xe1=0.05:0.05:4     %Plane 1
    for ze1=1:0.05:3          %ye1=0
      m_4=m_4+1;
   
   R91=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
   R93=sqrt((xe4-xe1)^2+(ye14-ye1)^2+(ze14-ze1)^2);
   R94=sqrt((xr-xe1)^2+(yr-ye1)^2+(zr-ze1)^2);
   R95=R91+R92+R93+R94;
  
  %angel for Second Ref
  theta4_11=acos(abs(zt-ze14)/R91);
  theta4_12=acos(2/R91);             %xe4-xt=4-2=2
  theta4_13=acos(abs(ze-ze1)/R93);            
  theta4_14=acos(abs(ze-ze14)/R93);
  theta4_15=acos(1/R94);             %yr-ye1=1-0=1
  theta4_16=acos(abs(zr-ze1)/R94);          
  
  Pres17=(((n+1)*(n+1))^2/(8*(pi)^3*(R91)^2*(R93)^2*(R94)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta4_11))^n*(cos(theta4_12))^n*cos(theta4_13)*(cos(theta4_14))^n*cos(theta4_15)*cos(theta4_16);
   T17=R95/S;
   Y_2(m_4)=Pres17;
   Z_2(m_4)=T17;
    end
       end
    end  
end

for ye14=0.05:0.05:8    %Plane 4
    for ze14=1:0.05:3    %xe4=4
   m_5=m_5+1;
   
    R96=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
    R97=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);

    %angel for First Ref
    theta4_17=acos(2/R96);
    phi37=acos(abs(zt-ze14)/R96);
    theta4_18=acos(3/R97);
    phi38=acos(abs(zr-ze14)/R97);
    
    Pre18=((n+1)*(n+1)/(4*(pi)^2*(R96)^2*(R97)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta4_17))^n*cos(phi37)*(cos(theta4_18))^n*cos(phi38);
    X_3(m_5)=Pre18;
  

    for xe12=0.05:0.05:4     %Plane 2
     for ze12=1:0.05:3        %ye2=8                                                                                             
         m_6=m_6+1;
   
   R96=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
   R98=sqrt((xe4-xe12)^2+(ye14-ye2)^2+(ze14-ze12)^2);
   R99=sqrt((xr-xe12)^2+(yr-ye2)^2+(zr-ze12)^2);
    R100=R96+R97+R98+R99;
    
    %angel for Second Ref
   theta4_19=acos(abs(zt-ze12)/R96);
   theta4_20=acos(2/R96);             %xe4-xt=4-2=2
   theta4_21=acos(abs(ze-ze14)/R98);     
   theta4_22=acos(abs(ze-ze12)/R98);
   theta4_23=acos(7/R99);             %yr-ye2=8-1=7
   theta4_24=acos(abs(zr-ze12)/R99);   
  
     
   Pres18=(((n+1)*(n+1))^2/(8*(pi)^3*(R96)^2*(R98)^2*(R99)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta4_19))^n*(cos(theta4_20))^n*cos(theta4_21)*(cos(theta4_22))^n*cos(theta4_23)*cos(theta4_24);
    T18=R100/S;
    Y_3(m_6)=Pres18;
    Z_3(m_6)=T18;
    
    end
       end
    end
   
end


   for ye14=0.05:0.05:8    %Plane 4
    for ze14=1:0.05:3    %xe4=4
   m_7=m_7+1;
   
    R101=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
    R102=sqrt((xr-xe4)^2+(yr-ye14)^2+(zr-ze14)^2);

    %angel for First Ref
    theta4_25=acos(2/R101);        %xe4-xt=4-2=2
    phi39=acos(abs(zt-ze14)/R101);
    theta4_26=acos(3/R102);        %xe4-xr=4-1=3
    phi40=acos(abs(zr-ze14)/R102);
    
    Pre19=((n+1)*(n+1)/(4*(pi)^2*(R101)^2*(R102)^2))*Tf*Tc*Pt*dA*Ar*(cos(theta4_25))^n*cos(phi39)*(cos(theta4_26))^n*cos(phi40);
    X_4(m_7)=Pre19;
  

    for ye13=0.05:0.05:8       %Plane 3
     for ze13=1:0.05:3          %xe3=0
       m_8=m_8+1;

     R101=sqrt((xt-xe4)^2+(yt-ye14)^2+(zt-ze14)^2);
     R103=sqrt((xe4-xe3)^2+(ye14-ye13)^2+(ze14-ze13)^2);
     R104=sqrt((xr-xe3)^2+(yr-ye13)^2+(zr-ze13)^2);
     R105=R101+R102+R103+R104;
      
      %angel for Second Ref
      theta4_27=acos(abs(zt-ze14)/R101);
      theta4_28=acos(2/R101);        %xe4-xt=4-2=2
      theta4_29=acos(abs(ze-ze14)/R103);        
      theta4_30=acos(abs(ze-ze13)/R103);
      theta4_31=acos(abs(zr-ze13)/R104);      
      theta4_32=acos(1/R104);         %xr-xe3=1-0=1 
      
    Pres19=(((n+1)*(n+1))^2/(8*(pi)^3*(R101)^2*(R103)^2*(R104)^2))*Tf*Tc*Pt*dA*dA2*Ar*(cos(theta4_27))^n*(cos(theta4_28))^n*cos(theta4_29)*(cos(theta4_30))^n*cos(theta4_31)*cos(theta4_32);
     T19=R105/S;
     Y_4(m_8)=Pres19; 
     Z_4(m_8)=T19;
 
     end
    end
   end
end
