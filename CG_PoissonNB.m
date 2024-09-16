function [f,B,S,se] = CG_PoissonNB(events,exposure,X,w,R,model,Beta)

if isequal(model,'Poisson')
    B = [log(sum(events.*w)./sum(exposure.*w))';zeros(size(X{1},2),1)];
else
    B = [log(sum(events.*w)./sum(exposure.*w))';zeros(size(X{1},2),1);0];
end

B(~isnan(Beta)) = Beta(~isnan(Beta));
sEL             = find(~isinf(B) & isnan(Beta));
delta           = 10^-5;
dL              = [];
for r = 1:R
    f(r) = LikeLiHood(B,events,exposure,X,w,model);
    if numel(B) ~= numel(dL)
        dL      = d_LikeLiHood(f(r),B,events,exposure,X,w,model);
    else
        dL(:,1) = dL(:,end);
        dL(:,2) = d_LikeLiHood(f(r),B,events,exposure,X,w,model);
    end
    
    if size(dL,2) > 1
        %b(1)    = dL(sEL,2)'*dL(sEL,2)*pinv(dL(sEL,1)'*dL(sEL,1));                      % Fletcher-Reeves
        %b(2)    = dL(sEL,2)'*(dL(sEL,1) - dL(sEL,2))*pinv(dL(sEL,1)'*dL(sEL,1));        % Polak-Ribière
        %b(3)    = dL(sEL,2)'*(dL(sEL,1) - dL(sEL,2))*pinv(-s'*(dL(sEL,1) - dL(sEL,2))); % Hestenes-Stiefel
        b(4)     = dL(sEL,2)'*dL(sEL,2)*pinv(-s'*(dL(sEL,1) - dL(sEL,2)));               % Dai-Yuan
        
        s        = -dL(sEL,2) + max(0,b(4))*s;
        for i = 1:2
            Bi      = B;
            Bi(sEL) = B(sEL) + delta*((-1)^i)*s;
            fi(i)   = LikeLiHood(Bi,events,exposure,X,w,model);
        end
        a       = -delta/2*(fi(2) - fi(1))*pinv((fi(2) - 2*f(end) + fi(1)));
        clear fi b i Bi
    else
        s       = -dL(sEL,end);
        a       = .01;
    end
    dL   = dL(:,end);
    
    for i = 1:10
        Bi      = B;
        Bi(sEL) = B(sEL) + a*s;
        fi(i)   = LikeLiHood(Bi,events,exposure,X,w,model);
        if fi(i) < f(r)
            a     = 0.75*a*(-1)^i;
        end
    end
    
    if (abs(a) < 10^-5 || (a*s)'*(a*s) < 10^-5 || abs(f(r) - fi(end)) < 10^-5) && r > 5
        break
    else
        clc;
        B(sEL)  = B(sEL) + a*s
    end
    clear Xi fi
end

H        = Hessian(B,events,exposure,X,w,model);
s        = find(isinf(B) | ~isnan(Beta));
H(:,s)   = NaN;
H(s,:)   = NaN;
s        = find(~isnan(H));
H        = reshape(H(s),numel(sEL),numel(sEL));
S        = pinv(-H);
se       = NaN(size(B));
se(sEL)  = sqrt(diag(S));


function D = d_LikeLiHood(f,B,events,exposure,X,w,model)
d = 10^-5;
D = NaN(size(B));
for i = 1:numel(B)
    Bi    = B;
    Bi(i) = Bi(i) + d;
    fi    = LikeLiHood(Bi,events,exposure,X,w,model);
    D(i)  = (fi - f)/d;
end

function lnL = LikeLiHood(B,events,exposure,X,w,model)
k = size(events,2);
for j = 1:size(B,2)
    if isequal(model(1:3),'NB1')
        alpha    = exp(B(end,j))';
        Beta     = B(1:end - 1,j);
        AeXB     = alpha.*exp(Beta(1:k)' + cell2mat(X)*kron(eye(k),Beta(k + 1:end))).*exposure(:,1:k);
        C        = gamma(events(:,1:k) + 1./alpha)./(gamma(events(:,1:k) + 1).*gamma(1./alpha));
        p        = C.*((1./(1 + AeXB)).^(1./alpha)).*((AeXB./(1 + AeXB)).^events(:,1:k));
        lnL(j,1) = sum(sum(w(:,1:k).*log(p)));
    elseif isequal(model,'NB2')
        delta    = exp(B(end,j))';
        Beta     = B(1:end - 1,j);
        mu       = exp(Beta(1:k)' + cell2mat(X)*kron(eye(k),Beta(k + 1:end))).*exposure(:,1:k);
        p        = (mu.^events(:,1:k)).*((1/(1 + delta)).^(mu/delta + events(:,1:k)));
        lnL(j,1) = sum(sum(w(:,1:k).*log(p)));
    else
        eXB      = exp(B(1:k,j)' + cell2mat(X)*kron(eye(k),B(k + 1:end,j))).*exposure(:,1:k);
        p        = (1./factorial(events(:,1:k))).*(eXB.^events(:,1:k)).*exp(-eXB);
        lnL(j,1) = sum(sum(w(:,1:k).*log(p)));
    end
end

function D = dlnLikelihood(B,events,exposure,X,w,model)
d = 10^-5;
for j = 1:size(B,2)
    for i = 1:size(B,1)
        Bi      = [B(:,j) B(:,j)];
        Bi(i,:) = Bi(i,:) + [d,-d];
        lnL     = LikeLiHood(Bi,events,exposure,X,w,model);
        D(i,j)  = (lnL(1) - lnL(2))/(2*d);
    end
end

function H = Hessian(B,events,exposure,X,w,model)
d = 10^-5;
for i = 1:numel(B)
    Bi      = [B B];
    Bi(i,:) = Bi(i,:) + [d,-d];
    dlnL    = dlnLikelihood(Bi,events,exposure,X,w,model);
    H(:,i)  = (dlnL(:,1) - dlnL(:,2))/(2*d);
end