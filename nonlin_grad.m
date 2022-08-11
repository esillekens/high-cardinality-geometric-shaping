function [scale,jvp] = nonlin_grad(c,x)

[M,D] = size(x);

% assert(D==2,'Only implemetend for 2D');

v = mean(sum(abs(x).^2,2).^2,1);
w = mean(sum(abs(x).^2,2),1).^2;

eK = v/w-2;

scale = ( 1 + c*eK ).^(-1/6);

% vprime = 4/M*x.*sum(abs(x).^2,2);
% wprime = 4/M*x.*mean(sum(abs(x).^2,2));


vwprime = 4/M*x.*(sum(abs(x).^2,2).*mean(sum(abs(x).^2,2)) ...
                    - mean(sum(abs(x).^2,2).^2))...
            ./mean(sum(abs(x).^2,2)).^3;
            
scale_grad = -1/6*scale.^7*c*vwprime;

jac = eye(M*D)*scale + x(:) * scale_grad(:).';

% grad_out = grad(:).'*jac;
% grad =reshape(grad_out,M,D);

jvp = @(grad) reshape( grad(:).'*jac , M,D);

end