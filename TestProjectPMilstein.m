 %Strong Convergence in terms of h
clc;
M=[2^8,2^9,2^10,2^11,2^12,2^13,2^14];
N=10^4;
T=1;
X0=1;
rep=10;

relerrorStrongFinal =zeros(rep,7);
relerrorStrongFinalT =zeros(rep,7);
relerrorStrongFinalT12 =zeros(rep,7);
relerrorStrongFinalRef =zeros(rep,7);


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

[X1FT,X1CT,X11FT,X11CT,X111FT,X111CT]=AdaptiveTamedEulerMilstein(M(1),N,T,8);
[X2FT,X2CT,X22FT,X22CT,X222FT,X222CT]=AdaptiveTamedEulerMilstein(M(2),N,T,9);
[X3FT,X3CT,X33FT,X33CT,X333FT,X333CT]=AdaptiveTamedEulerMilstein(M(3),N,T,10);
[X4FT,X4CT,X44FT,X44CT,X444FT,X444CT]=AdaptiveTamedEulerMilstein(M(4),N,T,11);
[X5FT,X5CT,X55FT,X55CT,X555FT,X555CT]=AdaptiveTamedEulerMilstein(M(5),N,T,12);
[X6FT,X6CT,X66FT,X66CT,X666FT,X666CT]=AdaptiveTamedEulerMilstein(M(6),N,T,13);
[X7FT,X7CT,X77FT,X77CT,X777FT,X777CT]=AdaptiveTamedEulerMilstein(M(7),N,T,14);

relerrorStrong = zeros(1,7);
relerrorStrongT = zeros(1,7);
relerrorStrongT12 = zeros(1,7);
relerrorStrongRef = zeros(1,7);


relerrorStrong2 = zeros(1,7);
relerrorStrongT2 = zeros(1,7);
relerrorStrongT122 = zeros(1,7);
relerrorStrongRef2 = zeros(1,7);

relerrorStrong3 = zeros(1,7);
relerrorStrongT3 = zeros(1,7);
relerrorStrongT133 = zeros(1,7);
relerrorStrongRef3 = zeros(1,7);


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

relerrorStrongT3(1)=relerrorStrongT3(1) + mean((X111FT(M(1)+1,:)-X111CT(M(1)./2+1,:)).^2 ); 
relerrorStrongT3(2)=relerrorStrongT3(2) + mean((X222FT(M(2)+1,:)-X222CT(M(2)./2+1,:)).^2 ); 
relerrorStrongT3(3)=relerrorStrongT3(3) + mean((X333FT(M(3)+1,:)-X333CT(M(3)./2+1,:)).^2 ); 
relerrorStrongT3(4)=relerrorStrongT3(4) + mean((X444FT(M(4)+1,:)-X444CT(M(4)./2+1,:)).^2 ); 
relerrorStrongT3(5)=relerrorStrongT3(5) + mean((X555FT(M(5)+1,:)-X555CT(M(5)./2+1,:)).^2 ); 
relerrorStrongT3(6)=relerrorStrongT3(6) + mean((X666FT(M(6)+1,:)-X666CT(M(6)./2+1,:)).^2 ); 
relerrorStrongT3(7)=relerrorStrongT3(7) + mean((X777FT(M(7)+1,:)-X777CT(M(7)./2+1,:)).^2 );


relerrorStrongFinal(i,:) = relerrorStrong;
relerrorStrongFinalT(i,:) = relerrorStrongT;
relerrorStrongFinalT12(i,:) = relerrorStrongT12;

relerrorStrongFinal2(i,:) = relerrorStrong2;
relerrorStrongFinalT2(i,:) = relerrorStrongT2;
relerrorStrongFinalT122(i,:) = relerrorStrongT122;

relerrorStrongFinal3(i,:) = relerrorStrong3;
relerrorStrongFinalT3(i,:) = relerrorStrongT3;
relerrorStrongFinalT133(i,:) = relerrorStrongT133;

%relerrorStrongFinalRef(i,:) = relerrorStrongRef;
end

relerrorFinalTrue = zeros(1,7);

%Tamed

relerrorFinalTrueT = zeros(1,7);

relerrorFinalTrue2 = zeros(1,7);

%Tamed

relerrorFinalTrueT2 = zeros(1,7);


relerrorFinalTrueT3 = zeros(1,7);

for i=1:7
  relerrorFinalTrueT2(i) = sqrt(mean(relerrorStrongFinalT2(:,i)));
end



plot(log2(M),log2(relerrorFinalTrueT2(1:7)),'-d')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueT(i) = sqrt(mean(relerrorStrongFinalT(:,i)));
end


plot(log2(M),log2(relerrorFinalTrueT(1:7)),'-x')

grid on;
hold on;



for i=1:7
  relerrorFinalTrueT3(i) = sqrt(mean(relerrorStrongFinalT3(:,i)));
end


plot(log2(M),log2(relerrorFinalTrueT3(1:7)),'-*')

grid on;
hold on;

plot(log2(M),log2(4./(M.^(1))),'--')

grid on;
hold on;

%legend('Adaptive','Tamed','slope -1', 'slope -0.5','location','northeast')
legend('Tamed Euler scheme (Example 1)', 'Tamed Milstein scheme (1) (Example 2)', 'Tamed Milstein scheme (2) (Example 2)', 'slope -1','location','northeast')

xlabel('Level l')
ylabel('log_2(RMSE)')

grid on;

%mean(mean(vecindices))
%mean(mean(vecindicesC))
%2^4

%mean(mean(vecindices2))
%mean(mean(vecindicesC2))
%2^5 


 

