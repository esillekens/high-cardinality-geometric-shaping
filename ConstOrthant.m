function x = ConstOrthant(x)

[M,D] = size(x);

bmap = de2bi((0:M-1).')==1;
orthant = ones(M,1,'logical');

[~,candidate] = max(abs(bmap.'*x));

flip =1-2*eye(D);
for d=1:D
 dist = x(bmap(:,candidate(d)),:)-flip(d,:).*x(~bmap(:,candidate(d)),:);
 if sum(abs(dist))<sum(M*eps(x))
     orthant = orthant&bmap(:,candidate(d));
 else
     warning('no symmetry')
 end
 
end

x = sign(mean(x(orthant,:))).*x(orthant,:);