%Fine Coarse Fine Coarse

function [X,X3,X1,X2,X4,X5] = AdaptiveTamedEulerMilstein(M,N,T,L)

  X0=1;

  h= T./M;
  X=zeros(M+1,N);
  X1=zeros(M+1,N);
  X4=zeros(M+1,N);
  
  X(1,:)=X0;
  X1(1,:)=X0;
  X4(1,:)=X0;

  sigma =1.5;
  c=0.5;
  
  [~,~,dW,~] = BrownianPathCorrL(T,L,0,0,0,N);
  
    for m=1:M
      mm=mean(X(m,:)./(1+X(m,:)));
      for n=1:N
          b= X(m,n)*sigma^2./2 -X(m,n)^3+c*mm;
          X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +(X(m,n)+mm)*dW(N+n,m) + 0.5*(X(m,n)+mm)*(dW(N+n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*dW(N+n,m) + 0.5*(X(m,n))*(dW(N+n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
    end
    
     for m=1:M
      mm2=mean(X1(m,:));
      for n=1:N
          b= X1(m,n)*sigma^2./2 -X1(m,n)^3+c*mm2;
          X1(m+1,n)= X1(m,n) + b*h./(1+ M^(-1)*abs(b)) + mm2*dW(N+n,m);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
     end
    
     for m=1:M
      mm3=mean(X4(m,:)./(1+X4(m,:)));
      for n=1:N
          b= X4(m,n)*sigma^2./2 -X4(m,n)^3+c*mm3;
          X4(m+1,n)= X4(m,n) + b*h./(1+ M^(-1)*abs(b)^2) +(X4(m,n)+mm3)*dW(N+n,m) + 0.5*(X4(m,n)+mm3)*(dW(N+n,m)^2-h);
          %X4(m+1,n)= X4(m,n) + b*h./(1+ M^(-1)*abs(b)^2) +X4(m,n)*dW(N+n,m) + 0.5*(X4(m,n))*(dW(N+n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
    end
     
  M1=M./2;
  h= T./M1;
  X3=zeros(M1+1,N);
  X3(1,:)=X0;
  
  X2=zeros(M1+1,N);
  X2(1,:)=X0; 
  
  X5=zeros(M1+1,N);
  X5(1,:)=X0; 
  
   for m=1:M1
      mmC=mean(X3(m,:)./(1+X3(m,:)));
      for n=1:N
          b= X3(m,n)*sigma^2./2 -X3(m,n)^3+c*mmC;
          X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) +(X3(m,n)+mmC)*dW(n,m) + 0.5*(X3(m,n)+mmC)*(dW(n,m)^2-h);
          %X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) +X3(m,n)*dW(n,m) + 0.5*(X3(m,n))*(dW(n,m)^2-h);
          %X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) +mmC*dW(n,m) ;
      end
   end
   
   for m=1:M1
      mmC3=mean(X5(m,:)./(1+X5(m,:)));
      for n=1:N
          b= X5(m,n)*sigma^2./2 -X5(m,n)^3+c*mmC3;
          X5(m+1,n)= X5(m,n) + b*h./(1+ M1^(-1)*abs(b)^2) +(X5(m,n)+mmC3)*dW(n,m) + 0.5*(X5(m,n)+mmC3)*(dW(n,m)^2-h);
          %X5(m+1,n)= X5(m,n) + b*h./(1+ M1^(-1)*abs(b)^2) +X5(m,n)*dW(n,m) + 0.5*(X5(m,n))*(dW(n,m)^2-h);
          %X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) +mmC*dW(n,m) ;
      end
   end
    
     
   for m=1:M1
      mmC2=mean(X2(m,:));
      for n=1:N
          b= X2(m,n)*sigma^2./2 -X2(m,n)^3+c*mmC2;
          X2(m+1,n)= X2(m,n) + b*h./(1+ M1^(-1)*abs(b)) + mmC2*dW(n,m);
          %X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) +mmC*dW(n,m) ;
      end
   end
       
end
