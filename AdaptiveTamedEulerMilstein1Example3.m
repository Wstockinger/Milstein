function [X,X3,X4,X5,X6,X7] = AdaptiveTamedEulerMilstein1Example3(M,N,T,L)

  X0=1;

  h= T./M;

  X=zeros(M+1,N);
  X4=zeros(M+1,N);
  X6=zeros(M+1,N);
  
  X(1,:)=X0;
  X6(1,:)=X0;
  X4(1,:)=X0;
  tmp1 = zeros(1,N);
  tmp2 = zeros(1,N);
  sigma =1.5;
  
  [~,~,dW,~] = BrownianPathCorrL(T,L,0,0,0,N);
  
  

    for m=1:M
        
      for n=1:N
        tmp1(n) = sum(sin(X(m,n)-X(m,:)));
        tmp2(n) = sum(cos(X(m,n)-X(m,:)));
      end
      for n=1:N
          b= X(m,n)*sigma^2./2 -X(m,n)^3+tmp1(n)./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) + (X(m,n)+tmp1(n)./N)*dW(N+n,m) + 0.5*(X(m,n)+tmp1(n)./N)*(1 + tmp2(n)./N)*(dW(N+n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
    end
    
    for m=1:M
      for n=1:N
          b= X6(m,n)*sigma^2./2 -X6(m,n)^3+sum(sin(X6(m,n)-X6(m,:)))./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X6(m+1,n)= X6(m,n) + b*h./(1+ M^(-1)*abs(b)) +(X6(m,n)+sum(sin(X6(m,n)-X6(m,:)))./N)*dW(N+n,m);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
    end
    
     for m=1:M
      for n=1:N
        tmp1(n) = sum(sin(X4(m,n)-X4(m,:)));
        tmp2(n) = sum(cos(X4(m,n)-X4(m,:)));
      end
      for n=1:N
          b= X4(m,n)*sigma^2./2 -X4(m,n)^3+tmp1(n)./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X4(m+1,n)= X4(m,n) + b*h./(1+ M^(-1)*abs(b)^2) + (X4(m,n)+tmp1(n)./N)*dW(N+n,m) + 0.5*(X4(m,n)+tmp1(n)./N)*(1 + tmp2(n)./N)*(dW(N+n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
     end
     
  M1=M./2;
  h= T./M1;
  X3=zeros(M1+1,N);
  X3(1,:)=X0;
  
  X7=zeros(M1+1,N);
  X7(1,:)=X0;
  
  X5=zeros(M1+1,N);
  X5(1,:)=X0; 

  
  
   for m=1:M1
       for n=1:N
        tmp1(n) = sum(sin(X3(m,n)-X3(m,:)));
        tmp2(n) = sum(cos(X3(m,n)-X3(m,:)));
      end
      for n=1:N
          b= X3(m,n)*sigma^2./2 -X3(m,n)^3+tmp1(n)./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X3(m+1,n)= X3(m,n) + b*h./(1+ M1^(-1)*abs(b)) + (X3(m,n)+tmp1(n)./N)*dW(n,m) + 0.5*(X3(m,n)+tmp1(n)./N)*(1 + tmp2(n)./N)*(dW(n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
   end
   
   for m=1:M1
      for n=1:N
        tmp1(n) = sum(sin(X5(m,n)-X5(m,:)));
        tmp2(n) = sum(cos(X5(m,n)-X5(m,:)));
      end
      for n=1:N
          b= X5(m,n)*sigma^2./2 -X5(m,n)^3+tmp1(n)./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X5(m+1,n)= X5(m,n) + b*h./(1+ M1^(-1)*abs(b)^2) + (X5(m,n)+tmp1(n)./N)*dW(n,m) + 0.5*(X5(m,n)+tmp1(n)./N)*(1 + tmp2(n)./N)*(dW(n,m)^2-h);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
   end
   
   for m=1:M1
      for n=1:N
          b= X7(m,n)*sigma^2./2 -X7(m,n)^3+sum(sin(X7(m,n)-X7(m,:)))./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X7(m+1,n)= X7(m,n) + b*h./(1+ M1^(-1)*abs(b)) +(X7(m,n)+sum(sin(X7(m,n)-X7(m,:)))./N)*dW(n,m);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
   end
   
     
end
