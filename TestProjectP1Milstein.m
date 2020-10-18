 %Strong Convergence in terms of h
clc;
M=[2^6,2^7,2^8,2^9,2^10,2^11,2^12];
%M=[2^3,2^4,2^5,2^6,2^7,2^8,2^9];
%N=20;
N=3;
T=1;
X0=1;
rep=500;

relerrorStrongFinal =zeros(rep,7);
relerrorStrongFinalT =zeros(rep,7);
relerrorStrongFinalT12 =zeros(rep,7);
relerrorStrongFinalRef =zeros(rep,7);
relerrorStrongFinalT4 = zeros(rep,7);

relerrorStrongFinal2 =zeros(rep,7);
relerrorStrongFinalT2 =zeros(rep,7);
relerrorStrongFinalT122 =zeros(rep,7);
relerrorStrongFinalRef2 =zeros(rep,7);
%[XRef] = AdaptiveTamedEuler2(14,N,T,14);

relerrorStrongFinal3 =zeros(rep,7);
relerrorStrongFinalT3 =zeros(rep,7);
relerrorStrongFinalT133 =zeros(rep,7);
relerrorStrongFinalRef3 =zeros(rep,7);

for i=1:rep

i

[X1FT,X1CT,X11FT,X11CT,YF,YC,ZF,ZC]=AdaptiveTamedEulerMilstein1(M(1),N,T,6);
[X2FT,X2CT,X22FT,X22CT,Y1F,Y1C,Z1F,Z1C]=AdaptiveTamedEulerMilstein1(M(2),N,T,7);
[X3FT,X3CT,X33FT,X33CT,Y2F,Y2C,Z2F,Z2C]=AdaptiveTamedEulerMilstein1(M(3),N,T,8);
[X4FT,X4CT,X44FT,X44CT,Y3F,Y3C,Z3F,Z3C]=AdaptiveTamedEulerMilstein1(M(4),N,T,9);
[X5FT,X5CT,X55FT,X55CT,Y4F,Y4C,Z4F,Z4C]=AdaptiveTamedEulerMilstein1(M(5),N,T,10);
[X6FT,X6CT,X66FT,X66CT,Y5F,Y5C,Z5F,Z5C]=AdaptiveTamedEulerMilstein1(M(6),N,T,11);
[X7FT,X7CT,X77FT,X77CT,Y6F,Y6C,Z6F,Z6C]=AdaptiveTamedEulerMilstein1(M(7),N,T,12);

relerrorStrong = zeros(1,7);
relerrorStrongT = zeros(1,7);
relerrorStrongT12 = zeros(1,7);
relerrorStrongRef = zeros(1,7);


relerrorStrong2 = zeros(1,7);
relerrorStrongT2 = zeros(1,7);
relerrorStrongT122 = zeros(1,7);
relerrorStrongRef2 = zeros(1,7);

relerrorStrongT3 = zeros(1,7);
relerrorStrongT4 = zeros(1,7);

%relerrorStrongT(1)=relerrorStrongT(1) + mean((X1FT(M(1)+1,:)-X1CT(M(1)./2+1,:)).^2) ; 
%relerrorStrongT(2)=relerrorStrongT(2) + mean((X2FT(M(2)+1,:)-X2CT(M(2)./2+1,:)).^2) ; 

relerrorStrongT(1)=relerrorStrongT(1) + mean((X1FT(M(1)+1,:)-X1CT(M(1)./2+1,:)).^2 ); 
relerrorStrongT(2)=relerrorStrongT(2) + mean((X2FT(M(2)+1,:)-X2CT(M(2)./2+1,:)).^2 ); 
relerrorStrongT(3)=relerrorStrongT(3) + mean((X3FT(M(3)+1,:)-X3CT(M(3)./2+1,:)).^2 ); 
relerrorStrongT(4)=relerrorStrongT(4) + mean((X4FT(M(4)+1,:)-X4CT(M(4)./2+1,:)).^2 ); 
relerrorStrongT(5)=relerrorStrongT(5) + mean((X5FT(M(5)+1,:)-X5CT(M(5)./2+1,:)).^2 ); 
relerrorStrongT(6)=relerrorStrongT(6) + mean((X6FT(M(6)+1,:)-X6CT(M(6)./2+1,:)).^2 ); 
relerrorStrongT(7)=relerrorStrongT(7) + mean((X7FT(M(7)+1,:)-X7CT(M(7)./2+1,:)).^2 );

relerrorStrongT2(1)=relerrorStrongT2(1) + mean((X11FT(M(1)+1,:)-X11CT(M(1)./2+1,:)).^2 ); 
relerrorStrongT2(2)=relerrorStrongT2(2) + mean((X22FT(M(2)+1,:)-X22CT(M(2)./2+1,:)).^2 ); 
relerrorStrongT2(3)=relerrorStrongT2(3) + mean((X33FT(M(3)+1,:)-X33CT(M(3)./2+1,:)).^2 ); 
relerrorStrongT2(4)=relerrorStrongT2(4) + mean((X44FT(M(4)+1,:)-X44CT(M(4)./2+1,:)).^2 ); 
relerrorStrongT2(5)=relerrorStrongT2(5) + mean((X55FT(M(5)+1,:)-X55CT(M(5)./2+1,:)).^2 ); 
relerrorStrongT2(6)=relerrorStrongT2(6) + mean((X66FT(M(6)+1,:)-X66CT(M(6)./2+1,:)).^2 ); 
relerrorStrongT2(7)=relerrorStrongT2(7) + mean((X77FT(M(7)+1,:)-X77CT(M(7)./2+1,:)).^2 );

relerrorStrongT3(1)=relerrorStrongT3(1) + mean((YF(M(1)+1,:)-YC(M(1)./2+1,:)).^2 ); 
relerrorStrongT3(2)=relerrorStrongT3(2) + mean((Y1F(M(2)+1,:)-Y1C(M(2)./2+1,:)).^2 ); 
relerrorStrongT3(3)=relerrorStrongT3(3) + mean((Y2F(M(3)+1,:)-Y2C(M(3)./2+1,:)).^2 ); 
relerrorStrongT3(4)=relerrorStrongT3(4) + mean((Y3F(M(4)+1,:)-Y3C(M(4)./2+1,:)).^2 ); 
relerrorStrongT3(5)=relerrorStrongT3(5) + mean((Y4F(M(5)+1,:)-Y4C(M(5)./2+1,:)).^2 ); 
relerrorStrongT3(6)=relerrorStrongT3(6) + mean((Y5F(M(6)+1,:)-Y5C(M(6)./2+1,:)).^2 ); 
relerrorStrongT3(7)=relerrorStrongT3(7) + mean((Y6F(M(7)+1,:)-Y6C(M(7)./2+1,:)).^2 );

relerrorStrongT4(1)=relerrorStrongT4(1) + mean((ZF(M(1)+1,:)-ZC(M(1)./2+1,:)).^2 ); 
relerrorStrongT4(2)=relerrorStrongT4(2) + mean((Z1F(M(2)+1,:)-Z1C(M(2)./2+1,:)).^2 ); 
relerrorStrongT4(3)=relerrorStrongT4(3) + mean((Z2F(M(3)+1,:)-Z2C(M(3)./2+1,:)).^2 ); 
relerrorStrongT4(4)=relerrorStrongT4(4) + mean((Z3F(M(4)+1,:)-Z3C(M(4)./2+1,:)).^2 ); 
relerrorStrongT4(5)=relerrorStrongT4(5) + mean((Z4F(M(5)+1,:)-Z4C(M(5)./2+1,:)).^2 ); 
relerrorStrongT4(6)=relerrorStrongT4(6) + mean((Z5F(M(6)+1,:)-Z5C(M(6)./2+1,:)).^2 ); 
relerrorStrongT4(7)=relerrorStrongT4(7) + mean((Z6F(M(7)+1,:)-Z6C(M(7)./2+1,:)).^2 );


relerrorStrongFinal(i,:) = relerrorStrong;
relerrorStrongFinalT(i,:) = relerrorStrongT;
relerrorStrongFinalT12(i,:) = relerrorStrongT12;

relerrorStrongFinal2(i,:) = relerrorStrong2;
relerrorStrongFinalT2(i,:) = relerrorStrongT2;
relerrorStrongFinalT122(i,:) = relerrorStrongT122;

relerrorStrongFinalT3(i,:) = relerrorStrongT3;
relerrorStrongFinalT4(i,:) = relerrorStrongT4;


%relerrorStrongFinalRef(i,:) = relerrorStrongRef;
end

relerrorFinalTrue = zeros(1,7);

%Tamed

relerrorFinalTrueT = zeros(1,7);

relerrorFinalTrue2 = zeros(1,7);

%Tamed

relerrorFinalTrueT2 = zeros(1,7);

%Euler

relerrorFinalTrueT3 = zeros(1,7);

relerrorFinalTrueT4 = zeros(1,7);

%for i=1:7
%  relerrorFinalTrueT3(i) = sqrt(mean(relerrorStrongFinalT3(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueT3(1:7)),'-d')

%grid on;
%hold on;


for i=1:7
  relerrorFinalTrueT(i) = sqrt(mean(relerrorStrongFinalT(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueT(1:7)),'-x')

grid on;
hold on;

%for i=1:7
%  relerrorFinalTrueT2(i) = sqrt(mean(relerrorStrongFinalT2(:,i)));
%end


%plot(log2(M),log2(relerrorFinalTrueT2(1:7)),'-*')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueT4(i) = sqrt(mean(relerrorStrongFinalT4(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueT4(1:7)),'-+')

grid on;
hold on;

plot(log2(M),log2(4./(M.^(1))),'--')

grid on;
hold on;

plot(log2(M),log2(1./(8*M.^(0.5))),':')

grid on;
hold on;


%legend('Adaptive','Tamed','slope -1', 'slope -0.5','location','northeast')
legend('Tamed Milstein scheme (1) (Example 3)', 'Full-Tamed Milstein scheme (1) (Example 3)', 'slope -1','slope -0.5','location','northeast')
%legend('Tamed Milstein scheme (1) (Example 3)', 'Full-Tamed Milstein scheme (1) (Example 3)', 'slope -1','location','northeast')
%legend('Tamed Euler scheme (Example 3)','Tamed Milstein scheme (1) (Example 3)', 'Tamed Milstein scheme (2) (Example 3)',  'Full-Tamed Milstein scheme (1) (Example 3)', 'slope -1','slope -0.5','location','northeast')

xlabel('Level l')
ylabel('log_2(RMSE)')

grid on;

%mean(mean(vecindices))
%mean(mean(vecindicesC))
%2^4

%mean(mean(vecindices2))
%mean(mean(vecindicesC2))
%2^5 


 

