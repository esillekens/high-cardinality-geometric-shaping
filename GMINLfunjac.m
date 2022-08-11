function [fk,gk] = GMINLfunjac(SNR,c,xk)
% orthant back to full
Morth = size(xk,1);
xk = OrthantConst(xk);

% change power for nonlinear fibre response
[nk,jvp] = nonlin_grad(c,xk);

% calculated GMI and grad
[fk,gk] = GMI_withgrad_nD(SNR+20*log10(nk),xk);

% chain rule for gradient
gk = jvp(gk);
gk = normalise_grad(xk,gk);

% extract orthant
gk = gk(1:Morth,:);
% flip signs
fk = -fk; gk=-gk; %maximise
end