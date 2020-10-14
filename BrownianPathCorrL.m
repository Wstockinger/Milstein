
function [W1,W2,dW1,dW2] = BrownianPathCorrL(T,L,S0,S02,rho,M)
  
    N=2^L;
    dW1 = zeros(2*M,N);
    W1 = zeros(2*M,N+1);
    W1(:,1) = S0;
    
    dW2 = zeros(2*M,N);
    W2 = zeros(2*M,N+1);
    W2(:,1) = S02;
    
    j=L+1;
    stepsize = T./(2^(j-1));
    for m=1:M
        for i = 2:(2^(j-1)+1)    
            jj=M+m;
            dW2(jj,i-1) = randn*sqrt(stepsize);
            dW1(jj,i-1) = rho*dW2(jj,i-1) + sqrt(1-rho^2)*randn*sqrt(stepsize);
            W1(jj,i) = W1(jj,i-1) + dW1(jj,i-1);
            W2(jj,i) = W2(jj,i-1) + dW2(jj,i-1);
        end
    end
    k=L;
    for m=1:M    
       for i = 2:(2^(k-1)+1) 
           jk=m;
           dW1(jk,i-1) = sum(dW1(M+m,((i-2)*2^(L-k+1)+1):((i-1)*2^(L-k+1))));
           dW2(jk,i-1) = sum(dW2(M+m,((i-2)*2^(L-k+1)+1):((i-1)*2^(L-k+1))));
           W1(jk,i)=  W1(jk,i-1) + dW1(jk,i-1);
           W2(jk,i)=  W2(jk,i-1) + dW2(jk,i-1);    
       end
    end 
end