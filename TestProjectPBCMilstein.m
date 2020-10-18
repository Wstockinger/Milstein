 %Strong Convergence in terms of h
clc;
M=[2^7,2^8,2^9,2^10,2^11,2^12,2^13];
N=5;
%N=5;
NN=[2,5,10];
T=1;
X0=1;
rep=40;

relerrorStrongFinal =zeros(rep,7);
relerrorStrongFinalT =zeros(rep,7);
relerrorStrongFinalTA =zeros(rep,7);
relerrorStrongFinalTB =zeros(rep,7);
relerrorStrongFinalTC =zeros(rep,7);
relerrorStrongFinalTD =zeros(rep,7);

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

[X1FT,X1CT,X11FT,X11CT]=TamedEulerMilstein1(M(1),NN(1),T,7);
[X2FT,X2CT,X22FT,X22CT]=TamedEulerMilstein1(M(2),NN(1),T,8);
[X3FT,X3CT,X33FT,X33CT]=TamedEulerMilstein1(M(3),NN(1),T,9);
[X4FT,X4CT,X44FT,X44CT]=TamedEulerMilstein1(M(4),NN(1),T,10);
[X5FT,X5CT,X55FT,X55CT]=TamedEulerMilstein1(M(5),NN(1),T,11);
[X6FT,X6CT,X66FT,X66CT]=TamedEulerMilstein1(M(6),NN(1),T,12);
[X7FT,X7CT,X77FT,X77CT]=TamedEulerMilstein1(M(7),NN(1),T,13);

[X1FT2,X1CT2,X11FT2,X11CT2]=TamedEulerMilstein1(M(1),NN(2),T,7);
[X2FT2,X2CT2,X22FT2,X22CT2]=TamedEulerMilstein1(M(2),NN(2),T,8);
[X3FT2,X3CT2,X33FT2,X33CT2]=TamedEulerMilstein1(M(3),NN(2),T,9);
[X4FT2,X4CT2,X44FT2,X44CT2]=TamedEulerMilstein1(M(4),NN(2),T,10);
[X5FT2,X5CT2,X55FT2,X55CT2]=TamedEulerMilstein1(M(5),NN(2),T,11);
[X6FT2,X6CT2,X66FT2,X66CT2]=TamedEulerMilstein1(M(6),NN(2),T,12);
[X7FT2,X7CT2,X77FT2,X77CT2]=TamedEulerMilstein1(M(7),NN(2),T,13);

[X1FT3,X1CT3,X11FT3,X11CT3]=TamedEulerMilstein1(M(1),NN(3),T,7);
[X2FT3,X2CT3,X22FT3,X22CT3]=TamedEulerMilstein1(M(2),NN(3),T,8);
[X3FT3,X3CT3,X33FT3,X33CT3]=TamedEulerMilstein1(M(3),NN(3),T,9);
[X4FT3,X4CT3,X44FT3,X44CT3]=TamedEulerMilstein1(M(4),NN(3),T,10);
[X5FT3,X5CT3,X55FT3,X55CT3]=TamedEulerMilstein1(M(5),NN(3),T,11);
[X6FT3,X6CT3,X66FT3,X66CT3]=TamedEulerMilstein1(M(6),NN(3),T,12);
[X7FT3,X7CT3,X77FT3,X77CT3]=TamedEulerMilstein1(M(7),NN(3),T,13);



relerrorStrong = zeros(1,7);
relerrorStrongT = zeros(1,7);
relerrorStrongRef = zeros(1,7);
relerrorStrong2 = zeros(1,7);
relerrorStrongT2 = zeros(1,7);
relerrorStrongT122 = zeros(1,7);
relerrorStrongRef2 = zeros(1,7);
relerrorStrongT3 = zeros(1,7);
relerrorStrongT4 = zeros(1,7);

relerrorStrongT12 = zeros(1,7);
relerrorStrongT22 = zeros(1,7);
relerrorStrongT13 = zeros(1,7);
relerrorStrongT23= zeros(1,7);
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

relerrorStrongT12(1)=relerrorStrongT12(1) + mean((X1FT2(M(1)+1,:)-X1CT2(M(1)./2+1,:)).^2 ); 
relerrorStrongT12(2)=relerrorStrongT12(2) + mean((X2FT2(M(2)+1,:)-X2CT2(M(2)./2+1,:)).^2 ); 
relerrorStrongT12(3)=relerrorStrongT12(3) + mean((X3FT2(M(3)+1,:)-X3CT2(M(3)./2+1,:)).^2 ); 
relerrorStrongT12(4)=relerrorStrongT12(4) + mean((X4FT2(M(4)+1,:)-X4CT2(M(4)./2+1,:)).^2 ); 
relerrorStrongT12(5)=relerrorStrongT12(5) + mean((X5FT2(M(5)+1,:)-X5CT2(M(5)./2+1,:)).^2 ); 
relerrorStrongT12(6)=relerrorStrongT12(6) + mean((X6FT2(M(6)+1,:)-X6CT2(M(6)./2+1,:)).^2 ); 
relerrorStrongT12(7)=relerrorStrongT12(7) + mean((X7FT2(M(7)+1,:)-X7CT2(M(7)./2+1,:)).^2 );

relerrorStrongT22(1)=relerrorStrongT22(1) + mean((X11FT2(M(1)+1,:)-X11CT2(M(1)./2+1,:)).^2 ); 
relerrorStrongT22(2)=relerrorStrongT22(2) + mean((X22FT2(M(2)+1,:)-X22CT2(M(2)./2+1,:)).^2 ); 
relerrorStrongT22(3)=relerrorStrongT22(3) + mean((X33FT2(M(3)+1,:)-X33CT2(M(3)./2+1,:)).^2 ); 
relerrorStrongT22(4)=relerrorStrongT22(4) + mean((X44FT2(M(4)+1,:)-X44CT2(M(4)./2+1,:)).^2 ); 
relerrorStrongT22(5)=relerrorStrongT22(5) + mean((X55FT2(M(5)+1,:)-X55CT2(M(5)./2+1,:)).^2 ); 
relerrorStrongT22(6)=relerrorStrongT22(6) + mean((X66FT2(M(6)+1,:)-X66CT2(M(6)./2+1,:)).^2 ); 
relerrorStrongT22(7)=relerrorStrongT22(7) + mean((X77FT2(M(7)+1,:)-X77CT2(M(7)./2+1,:)).^2 );

relerrorStrongT13(1)=relerrorStrongT13(1) + mean((X1FT3(M(1)+1,:)-X1CT3(M(1)./2+1,:)).^2 ); 
relerrorStrongT13(2)=relerrorStrongT13(2) + mean((X2FT3(M(2)+1,:)-X2CT3(M(2)./2+1,:)).^2 ); 
relerrorStrongT13(3)=relerrorStrongT13(3) + mean((X3FT3(M(3)+1,:)-X3CT3(M(3)./2+1,:)).^2 ); 
relerrorStrongT13(4)=relerrorStrongT13(4) + mean((X4FT3(M(4)+1,:)-X4CT3(M(4)./2+1,:)).^2 ); 
relerrorStrongT13(5)=relerrorStrongT13(5) + mean((X5FT3(M(5)+1,:)-X5CT3(M(5)./2+1,:)).^2 ); 
relerrorStrongT13(6)=relerrorStrongT13(6) + mean((X6FT3(M(6)+1,:)-X6CT3(M(6)./2+1,:)).^2 ); 
relerrorStrongT13(7)=relerrorStrongT13(7) + mean((X7FT3(M(7)+1,:)-X7CT3(M(7)./2+1,:)).^2 );

relerrorStrongT23(1)=relerrorStrongT23(1) + mean((X11FT3(M(1)+1,:)-X11CT3(M(1)./2+1,:)).^2 ); 
relerrorStrongT23(2)=relerrorStrongT23(2) + mean((X22FT3(M(2)+1,:)-X22CT3(M(2)./2+1,:)).^2 ); 
relerrorStrongT23(3)=relerrorStrongT23(3) + mean((X33FT3(M(3)+1,:)-X33CT3(M(3)./2+1,:)).^2 ); 
relerrorStrongT23(4)=relerrorStrongT23(4) + mean((X44FT3(M(4)+1,:)-X44CT3(M(4)./2+1,:)).^2 ); 
relerrorStrongT23(5)=relerrorStrongT23(5) + mean((X55FT3(M(5)+1,:)-X55CT3(M(5)./2+1,:)).^2 ); 
relerrorStrongT23(6)=relerrorStrongT23(6) + mean((X66FT3(M(6)+1,:)-X66CT3(M(6)./2+1,:)).^2 ); 
relerrorStrongT23(7)=relerrorStrongT23(7) + mean((X77FT3(M(7)+1,:)-X77CT3(M(7)./2+1,:)).^2 );



relerrorStrongFinalT(i,:) = relerrorStrongT;
relerrorStrongFinalT2(i,:) = relerrorStrongT2;

relerrorStrongFinalTA(i,:) = relerrorStrongT12;
relerrorStrongFinalTB(i,:) = relerrorStrongT22;

relerrorStrongFinalTC(i,:) = relerrorStrongT13;
relerrorStrongFinalTD(i,:) = relerrorStrongT23;



end

relerrorFinalTrue = zeros(1,7);

%Tamed

relerrorFinalTrueT = zeros(1,7);

relerrorFinalTrue2 = zeros(1,7);

%Tamed

relerrorFinalTrueT2 = zeros(1,7);

relerrorFinalTrueTA = zeros(1,7);
relerrorFinalTrueTB = zeros(1,7);
relerrorFinalTrueTC = zeros(1,7);
relerrorFinalTrueTD = zeros(1,7);

%Euler

%relerrorFinalTrueT3 = zeros(1,7);

%relerrorFinalTrueT4 = zeros(1,7);

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

for i=1:7
  relerrorFinalTrueT2(i) = sqrt(mean(relerrorStrongFinalT2(:,i)));
end


plot(log2(M),log2(relerrorFinalTrueT2(1:7)),'-*')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueT4(i) = sqrt(mean(relerrorStrongFinalT4(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueT4(1:7)),'-+')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueTA(i) = sqrt(mean(relerrorStrongFinalTA(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueTA(1:7)),'-+')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueTB(i) = sqrt(mean(relerrorStrongFinalTB(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueTB(1:7)),'-o')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueTC(i) = sqrt(mean(relerrorStrongFinalTC(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueTC(1:7)),'-s')

grid on;
hold on;

for i=1:7
  relerrorFinalTrueTD(i) = sqrt(mean(relerrorStrongFinalTD(:,i)));
end

plot(log2(M),log2(relerrorFinalTrueTD(1:7)),'-d')

grid on;
hold on;

plot(log2(M),log2(4./(M.^(1))),'--')

grid on;
hold on;

plot(log2(M),log2(1./(8*(M.^(0.5)))),':')

grid on;
hold on;


%legend('Adaptive','Tamed','slope -1', 'slope -0.5','location','northeast')
legend('Tamed Euler scheme (1) (Example 1), N=2', 'Full-Tamed Milstein scheme (1) (Example 1), N=2', 'Tamed Euler scheme (1) (Example 1), N=5', 'Full-Tamed Milstein scheme (1) (Example 1), N=5','Tamed Euler scheme (1) (Example 1), N=10', 'Full-Tamed Milstein scheme (1) (Example 1), N=10' , 'slope -1','slope -0.5','location','northeast')
%legend('Tamed Euler scheme (Example 3)','Tamed Milstein scheme (1) (Example 3)', 'Tamed Milstein scheme (2) (Example 3)',  'Full-Tamed Milstein scheme (1) (Example 3)', 'slope -1','slope -0.5','location','northeast')

xlabel('Level l')
ylabel('log_2(RMSE)')

%ylim([-2 -13])

grid on;



 

