
%Strong Convergence in terms of h
clc;
M=[2^4,2^5,2^6];
%N=5;
%N=5;
NN=[2^2,2^3,2^4,2^5,2^6];
T=1;
X0=1;
rep=100;

relerrorStrongFinal =zeros(rep,5);
relerrorStrongFinalT =zeros(rep,5);
relerrorStrongFinalTA =zeros(rep,5);
relerrorStrongFinalTB=zeros(rep,5);

relerrorStrongFinalT12 =zeros(rep,5);
relerrorStrongFinalRef =zeros(rep,5);


relerrorStrongFinal2 =zeros(rep,5);
relerrorStrongFinalT2 =zeros(rep,5);
relerrorStrongFinalT122 =zeros(rep,5);
relerrorStrongFinalRef2 =zeros(rep,5);
%[XRef] = AdaptiveTamedEuler2(14,N,T,14);

relerrorStrongFinal3 =zeros(rep,5);
relerrorStrongFinalT3 =zeros(rep,5);
relerrorStrongFinalT133 =zeros(rep,5);
relerrorStrongFinalRef3 =zeros(rep,5);

for i=1:rep

i

[X1FT,X1CT,X11FT,X11CT]=TamedEulerMilsteinRevisionPaper1(M(1),NN(1),T,4);
[X2FT,X2CT,X22FT,X22CT]=TamedEulerMilsteinRevisionPaper1(M(1),NN(2),T,4);
[X3FT,X3CT,X33FT,X33CT]=TamedEulerMilsteinRevisionPaper1(M(1),NN(3),T,4);
[X4FT,X4CT,X44FT,X44CT]=TamedEulerMilsteinRevisionPaper1(M(1),NN(4),T,4);
[X5FT,X5CT,X55FT,X55CT]=TamedEulerMilsteinRevisionPaper1(M(1),NN(5),T,4);

[X1FTA,X1CTA,X11FTA,X11CTA]=TamedEulerMilsteinRevisionPaper1(M(2),NN(1),T,5);
[X2FTA,X2CTA,X22FTA,X22CTA]=TamedEulerMilsteinRevisionPaper1(M(2),NN(2),T,5);
[X3FTA,X3CTA,X33FTA,X33CTA]=TamedEulerMilsteinRevisionPaper1(M(2),NN(3),T,5);
[X4FTA,X4CTA,X44FTA,X44CTA]=TamedEulerMilsteinRevisionPaper1(M(2),NN(4),T,5);
[X5FTA,X5CTA,X55FTA,X55CTA]=TamedEulerMilsteinRevisionPaper1(M(2),NN(5),T,5);

[X1FTB,X1CTB,X11FTB,X11CTB]=TamedEulerMilsteinRevisionPaper1(M(3),NN(1),T,6);
[X2FTB,X2CTB,X22FTB,X22CTB]=TamedEulerMilsteinRevisionPaper1(M(3),NN(2),T,6);
[X3FTB,X3CTB,X33FTB,X33CTB]=TamedEulerMilsteinRevisionPaper1(M(3),NN(3),T,6);
[X4FTB,X4CTB,X44FTB,X44CTB]=TamedEulerMilsteinRevisionPaper1(M(3),NN(4),T,6);
[X5FTB,X5CTB,X55FTB,X55CTB]=TamedEulerMilsteinRevisionPaper1(M(3),NN(5),T,6);




relerrorStrong = zeros(1,5);
relerrorStrongT = zeros(1,5);
relerrorStrongTA = zeros(1,5);
relerrorStrongTB = zeros(1,5);

relerrorStrongRef = zeros(1,5);
relerrorStrong2 = zeros(1,5);
relerrorStrongT2 = zeros(1,5);
relerrorStrongT122 = zeros(1,5);
relerrorStrongRef2 = zeros(1,5);
relerrorStrongT3 = zeros(1,5);
relerrorStrongT4 = zeros(1,5);

relerrorStrongT12 = zeros(1,5);
relerrorStrongT22 = zeros(1,5);
relerrorStrongT13 = zeros(1,5);
relerrorStrongT23= zeros(1,5);
%relerrorStrongT(1)=relerrorStrongT(1) + mean((X1FT(M(1)+1,:)-X1CT(M(1)./2+1,:)).^2) ; 
%relerrorStrongT(2)=relerrorStrongT(2) + mean((X2FT(M(2)+1,:)-X2CT(M(2)./2+1,:)).^2) ; 

relerrorStrongT(1)=relerrorStrongT(1) + mean((X1FT(M(1)+1,:)-X11FT(M(1)+1,:)).^2 ); 
relerrorStrongT(2)=relerrorStrongT(2) + mean((X2FT(M(1)+1,:)-X22FT(M(1)+1,:)).^2 ); 
relerrorStrongT(3)=relerrorStrongT(3) + mean((X3FT(M(1)+1,:)-X33FT(M(1)+1,:)).^2 ); 
relerrorStrongT(4)=relerrorStrongT(4) + mean((X4FT(M(1)+1,:)-X44FT(M(1)+1,:)).^2 ); 
relerrorStrongT(5)=relerrorStrongT(5) + mean((X5FT(M(1)+1,:)-X55FT(M(1)+1,:)).^2 ); 

relerrorStrongTA(1)=relerrorStrongTA(1) + mean((X1FTA(M(2)+1,:)-X11FTA(M(2)+1,:)).^2 ); 
relerrorStrongTA(2)=relerrorStrongTA(2) + mean((X2FTA(M(2)+1,:)-X22FTA(M(2)+1,:)).^2 ); 
relerrorStrongTA(3)=relerrorStrongTA(3) + mean((X3FTA(M(2)+1,:)-X33FTA(M(2)+1,:)).^2 ); 
relerrorStrongTA(4)=relerrorStrongTA(4) + mean((X4FTA(M(2)+1,:)-X44FTA(M(2)+1,:)).^2 ); 
relerrorStrongTA(5)=relerrorStrongTA(5) + mean((X5FTA(M(2)+1,:)-X55FTA(M(2)+1,:)).^2 ); 

relerrorStrongTB(1)=relerrorStrongTB(1) + mean((X1FTB(M(3)+1,:)-X11FTB(M(3)+1,:)).^2 ); 
relerrorStrongTB(2)=relerrorStrongTB(2) + mean((X2FTB(M(3)+1,:)-X22FTB(M(3)+1,:)).^2 ); 
relerrorStrongTB(3)=relerrorStrongTB(3) + mean((X3FTB(M(3)+1,:)-X33FTB(M(3)+1,:)).^2 ); 
relerrorStrongTB(4)=relerrorStrongTB(4) + mean((X4FTB(M(3)+1,:)-X44FTB(M(3)+1,:)).^2 ); 
relerrorStrongTB(5)=relerrorStrongTB(5) + mean((X5FTB(M(3)+1,:)-X55FTB(M(3)+1,:)).^2 ); 

relerrorStrongFinalT(i,:) = relerrorStrongT;
relerrorStrongFinalT2(i,:) = relerrorStrongT2;

relerrorStrongFinalTA(i,:) = relerrorStrongTA;
relerrorStrongFinalTB(i,:) = relerrorStrongTB;
end

relerrorFinalTrue = zeros(1,7);

%Tamed

relerrorFinalTrueT = zeros(1,5);

relerrorFinalTrue2 = zeros(1,5);

%Tamed

relerrorFinalTrueT2 = zeros(1,5);

relerrorFinalTrueTA = zeros(1,5);
relerrorFinalTrueTB = zeros(1,5);
relerrorFinalTrueTC = zeros(1,5);
relerrorFinalTrueTD = zeros(1,5);

%Euler

%relerrorFinalTrueT3 = zeros(1,7);

%relerrorFinalTrueT4 = zeros(1,7);

%for i=1:7
%  relerrorFinalTrueT3(i) = sqrt(mean(relerrorStrongFinalT3(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueT3(1:7)),'-d')

%grid on;
%hold on;


for i=1:5
  relerrorFinalTrueT(i) = sqrt(mean(relerrorStrongFinalT(:,i)));
end

plot(log2(NN),log2(relerrorFinalTrueT(1:5)),'-x')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueT2(i) = sqrt(mean(relerrorStrongFinalT2(:,i)));
%end


%plot(log2(M),log2(relerrorFinalTrueT2(1:7)),'-*')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueT4(i) = sqrt(mean(relerrorStrongFinalT4(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueT4(1:7)),'-+')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueTA(i) = sqrt(mean(relerrorStrongFinalTA(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueTA(1:7)),'-+')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueTB(i) = sqrt(mean(relerrorStrongFinalTB(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueTB(1:7)),'-o')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueTC(i) = sqrt(mean(relerrorStrongFinalTC(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueTC(1:7)),'-s')

%grid on;
%hold on;

%for i=1:7
%  relerrorFinalTrueTD(i) = sqrt(mean(relerrorStrongFinalTD(:,i)));
%end

%plot(log2(M),log2(relerrorFinalTrueTD(1:7)),'-d')

grid on;
hold on;

%plot(log2(M),log2(4./(M.^(1))),'--')

for i=1:5
  relerrorFinalTrueTA(i) = sqrt(mean(relerrorStrongFinalTA(:,i)));
end

plot(log2(NN),log2(relerrorFinalTrueTA(1:5)),'-s')

grid on;
hold on;

for i=1:5
  relerrorFinalTrueTB(i) = sqrt(mean(relerrorStrongFinalTB(:,i)));
end

plot(log2(NN),log2(relerrorFinalTrueTB(1:5)),'-d')

grid on;
hold on;


plot(log2(NN),log2(1./(32*(NN.^(0.5)))),':')

grid on;
hold on;


legend('Convergence of the L-derivative terms (M=16)', 'Convergence of the L-derivative terms (M=32)','Convergence of the L-derivative terms (M=64)','slope -0.5','location','northeast')
%legend('Tamed Euler Scheme (1) (Example 1), N=2', 'Full-Tamed Milstein scheme (1) (Example 1), N=2', 'Tamed Euler Scheme (1) (Example 1), N=5', 'Full-Tamed Milstein scheme (1) (Example 1), N=5','Tamed Euler Scheme (1) (Example 1), N=10', 'Full-Tamed Milstein scheme (1) (Example 1), N=10' , 'slope -1','slope -0.5','location','northeast')
%legend('Tamed Euler scheme (Example 3)','Tamed Milstein scheme (1) (Example 3)', 'Tamed Milstein scheme (2) (Example 3)',  'Full-Tamed Milstein scheme (1) (Example 3)', 'slope -1','slope -0.5','location','northeast')

xlabel('Level l')
ylabel('log_2(RMSE)')

grid on;

