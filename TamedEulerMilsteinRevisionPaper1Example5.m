%Fine Coarse Fine Coarse

function [X,X8] = TamedEulerMilsteinRevisionPaper1Example5(M,N,T,L)

  X0=1;

  h= T./M;
  htemp = sqrt(2./h);
  X=zeros(M+1,N);
  XX = zeros(M+1,N);
  X8=zeros(M+1,N);
  XX8 = zeros(M+1,N);
  X8(1,:)=X0;
  X(1,:)=X0;
  tmp1 = zeros(1,N);
  tmp2 = zeros(1,N);
  tmp3 = zeros(1,N);
  sigma =1.5;
  
  [~,~,dW,~] = BrownianPathCorrL(T,L,0,0,0,N);
  
  levyarea = zeros(N*N,M);
  for m=1:M
      for i=1:N 
          for j=1:N
              sumhelp=0;
               for ii=1:2^7
                      %stdd = delta/(2*pi^2*ii^2);
                      randnn=randn(1,4);
                      sumhelp = sumhelp + (1./ii)*(randnn(1)*(randnn(2) - htemp*dW(N+j,m)) - randnn(3)*(randnn(4)-htemp*dW(N+i,m)));
               end
             levyarea(j+(i-1)*N,m) = - (h./(2*pi))*sumhelp;     
          end
      end
  end
  
   for m=1:M
      mm=mean(X8(m,:)); 
      XX8(m,:) = X8(m,:)./(1+X8(m,:).^2); 
      meang = mean(XX8(m,:));
      tmpm= sum((XX8(m,:) - meang).^2)./N;
      for n=1:N
        tmp1(n) = - (X8(m,n).^2-1)./((1+X8(m,n).^2).^2);
      end
      tmpM = exp(-tmpm);
      temp=0;
      for n=1:N
          b= X8(m,n)*sigma^2./2 -X8(m,n)^3+mm;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          levy = 0;
          for l=1:N
           if(n==l)
             tt=h;  
             levytemp = 0;
           else
             tt=0;
             levytemp = levyarea((n+(n-1)*N -temp) + l-1,m);
           end
           levy =  levy + (2*tmp1(l)*meang-2*XX8(m,l)*tmp1(l))*tmpM*tmpM*(0.5*dW(N+n,m)*dW(N+l,m) -0.5*tt + levytemp);
          end 
          X8(m+1,n)= X8(m,n) + b*h./(1+ M^(-1)*abs(b)) + tmpM*dW(N+n,m) + levy./N;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
          temp = temp +1;
      end
    end
    
  
  
    for m=1:M
      mm=mean(X(m,:)); 
      XX(m,:) = X(m,:)./(1+X(m,:).^2); 
      meang = mean(XX(m,:));
      tmp=sum((XX(m,:) - meang).^2)./N;
      %for n=1:N
      %  tmp1(n) = sum(sin(X(m,n)+X(m,:)));
      %end
      tmpM = exp(-tmp);
      for n=1:N
          %mm = X(m,n);
          b= X(m,n)*sigma^2./2 -X(m,n)^3+mm;
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +X(m,n)*mm*dW(N+n,m) + 0.5*(X(m,n)*mm^2)*(dW(N+n,m)^2-h);
          X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) + tmpM*dW(N+n,m);
          %X(m+1,n)= X(m,n) + b*h./(1+ M^(-1)*abs(b)) +mm*dW(N+n,m) ;
      end
    end
     
  