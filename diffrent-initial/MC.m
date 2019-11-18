%======================
%In this script you can calculate the energy for different initial states using Balance theory
%=========================
% clear
N =50;
T =0:0.1:N;
iter=50;
m = zeros(length(T), 2);
Energy=zeros(iter,length(T));


for t = 1:length(T)
    l = l+1;
    if mod(l,100) == 0
        l;
    end;
    for repeat=1:iter
   
         W = triu(sign(randn(N, N)));
         W = W+W';
         W = W-diag(diag(W));
  
        
        w = conflict_MC(W(1:end), T(t));
        m(l, :) = [w(1) mean(w)];
        rr=m(l,1)';
        Energy(repeat,l)=-rr;
        
    end;
    EnergyMean=mean(Energy);
end;

plot(T,EnergyMean,'.-');
xlabel('Temprature');
ylabel('Energy');

