% clear
N =50;
T =0:0.1:N;
l = 0;
iter=50;
m = zeros(length(T), 2);
Energy=zeros(iter,length(T));


for t = 1:length(T)
    l = l+1;
    if mod(l,100) == 0
        l;
    end;
    for repeat=1:iter
%         
        %W=-(ones(N)-eye(N));
%         W(1:N/2,1:N/2)=1;
%         W(N/2+1:N,N/2+1:N)=1;
         W = triu(sign(randn(N, N)));
         W = W+W';
         W = W-diag(diag(W));
   
%         
       % nNeg=5*repeat;
       % cc=find(triu(W==1));
       % inds=randsample(cc,nNeg);
       % W(inds)=-1;
       % W(W'==-1)=-1;
        % %         %======================
      %  nPos=5*repeat;
      % dd=find(triu(W==-1));
       %inds=randsample(dd,nPos);
       % W(inds)=1;
       % W(W'==1)=1;
        
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

% hRng = -1:0.05:0;
% dd = zeros(length(T), length(hRng));
% for jj=1:length(T)
%  
%  dd(jj,:)=hist(Energy(:,jj), hRng);
% end;
%  imagesc(T, -1:0.01:0, log10(dd)')
% surf(T, hRng, (dd)')
%  mesh(T, hRng, (dd)')
% contour(T, hRng, (dd)')   final graph is contour 