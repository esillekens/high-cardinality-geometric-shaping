function grad = normalise_grad(x,grad)

[M,D] = size(grad);

% all dimensions are identical
x = x(:);
grad = grad(:);

% calculate jacobian matrix and apply chain rule
% grad = (grad'*((x'*x*eye(D*M)-x*x')/(sqrt(1/M)*(x'*x)^1.5)))';
grad = ((grad'*(x'*x) - (grad'*x)*x')/(sqrt(1/M)*(x'*x)^1.5))';

% reshape to D dimensions
grad = reshape(grad,M,D);

end