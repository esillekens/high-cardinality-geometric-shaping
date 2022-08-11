function x = OrhtantConst(x)

[~,D] = size(x);

flip = 1-2*eye(D);
for d=1:D
    x = [x;flip(d,:).*x];
end