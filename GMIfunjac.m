function [fk,gk] = GMIfunjac(SNR,xk)
% number of points in orthant
Morth = size(xk,1);
% Orthant back to full constellation
xk = OrthantConst(xk);
% function and grad from GMI
[fk,gk] = GMI_withgrad_nD(SNR,xk);
% chain rule for normalisation
gk = normalise_grad(xk,gk);
% extract orthant
gk = gk(1:Morth,:);
% flip signs
fk = -fk; gk=-gk; %maximise
end