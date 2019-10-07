%clear all;
N =50;
% link=N*(N-1)/2;
Repeats=100;
numberOfIteration=50000;
T =1:1:40;
% T=10:1:26;
l = 0;
level0=zeros(Repeats,2);
bk0 = zeros(Repeats, length(T));
bk = zeros(Repeats, length(T));

p = zeros(Repeats*length(T), 3);
% E = zeros(length(T), 3);
for t =1:length(T)
    temperature=T(t);
    l = l+1;
    if mod(l, 100)== 0
        l;
    end;
    jj=0;
    ii=0;
    for iter = 1:Repeats
      
%          W = triu(sign(randn(N, N)));
%          nNeg=randi(500);
%          cc=find(triu(W==-1));
%        inds=randsample(cc,nNeg);
%        W(inds)=1;
%         W(W'==1)=1;
%        W = W-diag(diag(W));
        r = rand^2;
        W = rand(N)>r;
        W = sign(W-.5);
        W = triu(W, 1)+triu(W, 1)';
      
        bk0(iter,t) = (trace(W^3))./(6*nchoosek(N, 3));
%         bk0(iter, t) = 
        
        [w, e] = balance_MC(W(1:end), temperature, numberOfIteration);
        w = reshape(w, N, N);
        bk(iter,t) = (trace(w^3))./(6*nchoosek(N, 3));
        
       
    
        p(iter+(t-1)*Repeats, :) = [temperature -bk0(iter, t) -bk(iter, t)];
    end
    
%     E1 = floor((-bk0+1)*100);
%     v = zeros(200, 1);
%     for i = 1:100
%         v(i) = var(bk(E1==i));
%     end
%     [m, ix] = max(v);
%     
%     E(t, :) = [-mean(bk(E1>150)) -ix/100 -mean(bk(E1<10))];
%     
%     plot(-1:0.01:-0.01, v(1:100))
%     getframe;
% %     pause(.1)
   
    
end

% plot(T, E, 'o-')
% scatter(-bk0,-bk)
scatter3(p(:, 1), p(:, 2), p(:, 3))
xlabel('T')
ylabel('Initail <s s s>');
zlabel('Final <s s s>');

    
% b = glmfit(-bk0, [bk ones(length(bk), 1)], 'binomial', 'link', 'logit');
% yFit = glmval(b, -1:.01:1, 'logit')
% hold on
% plot(-1:.01:1, -yFit, 'r')




