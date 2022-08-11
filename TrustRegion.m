function [xopt,fhist] = TrustRegion(funjac, x0, niter, min_radius, verbose)

xk = x0;
Bk = eye(numel(xk));
Dk = 1;

% small numbers
eta = 1e-3; % \in (0,1e-3)
r = 0.3; % \in (0,1)


if nargin<3
niter = 100;
end
if nargin<4
    min_radius = 1e-6;
end
if nargin<5
    verbose = 1;
end

fhist = zeros(niter,1);

[fk,gk] = funjac(xk);

if verbose>=1
    fprintf('|%10s','Function')
    fprintf('|%11s','Act. Red')
    fprintf('|%10s','Act/Pred')
    fprintf('|%10s','Grad. mag')
    fprintf('|%10s','Trust.R')
    fprintf('|%10s','Step mag')
    fprintf('|\n')
end

kopt = 0;
xopt = xk;

for k = 1:niter
    
sk = tr_step(gk(:),Bk,Dk);

sk = sk + (abs(xk(:))-xk(:));

[fkp1,gkp1] = funjac( xk + reshape(sk,size(xk)) );

if all(fkp1<fhist)
    kopt = k;
    xopt = xk + reshape(sk,size(xk));
end

yk = gkp1(:)-gk(:);

ared = fk-fkp1; % actual reduction
pred = -( gkp1(:).'*sk + 1/2* sk.'*Bk*sk); %predicted reduction

if verbose>=1
    fprintf('|%10s',num2str(fk));
    fprintf('|%11s',num2str(ared));
    fprintf('|%10s',num2str(ared/pred));
    fprintf('|%10s',num2str(sqrt(sum(gk.^2,'all'))));
    fprintf('|%10s',num2str(Dk));
    fprintf('|%10s',num2str(sqrt(sum(abs(sk).^2,'all'))));
    fprintf('|\n')
end
fhist(k) = fk;

% if the step is good
if (ared/pred > eta) || (ared>0) % only update if the actual and prediction align
    xk = xk+reshape(sk,size(xk));
    fk = fkp1;
    gk = gkp1;
end

% adjust step size
step = (1+sqrt(5))/2;
if ared/pred > 0.8
    if sqrt(sum(abs(sk).^2,'all')) <= 0.8*Dk
        Dk = Dk;
    else
        Dk = 3*step*Dk;
    end
elseif 0.2<=ared/pred
    Dk = Dk;
else
    Dk = 1/step*Dk;
end

if Dk>numel(sk)
    Dk = numel(sk);
end

% only update if denominator is big enough
if abs(sk.'*(yk-Bk*sk)) >= r*sum(abs(sk).^2)*sum(abs(yk-Bk*sk).^2)
%     disp(abs(sk.'*(yk-Bk*sk))./(sum(abs(sk).^2)*sum(abs(yk-Bk*sk).^2)))
    Bk = Bk+ (yk-Bk*sk)*(yk-Bk*sk).'/((yk-Bk*sk).'*sk);
    % update even if xkp1 = xk
else
    Bk = Bk;
end

if Dk <= min_radius
if verbose>=1
    disp('Stopping reached radius')
end
    break
end

end
if verbose>=1
fprintf('Best iteration %d in %d iterations\n', kopt, k)
end