close all
clear

M = 256;
D = 2;

SNR = 15;

% initilise with a single orthant
X_init = randn(M/(2^D),D);
% relabel 
X_init = RelabelAPSK_spectral(X_init);

%% prepare optimisation
funjac = @(x) GMIfunjac(SNR, x);

%% plot initial constellation
[GMI_init, grad_init] = funjac(X_init);
gradient_quiver(X_init,grad_init)
title(['Initial Constellation GMI:', num2str(-GMI_init)])

%% Optimisation loop

max_iter = 1000;
min_region = 1e-6;
verbose = 1;

[X_opt, fun_step] = TrustRegion(funjac, X_init, max_iter, min_region, verbose);

%% plot result

[GMI_opt, grad_opt] = funjac(X_opt);
gradient_quiver(X_opt,grad_opt)
title(['Optimised Constellation GMI:', num2str(-GMI_opt)])


%% plot function
function gradient_quiver(X,G)

%convert Orthant to full
X = OrthantConst(X);
G = OrthantConst(G);

figure,
hold on
quiver(X(:,1),X(:,2),G(:,1),G(:,2))
plot(X(:,1),X(:,2),'.','MarkerSize',6)
grid on
axis equal

end