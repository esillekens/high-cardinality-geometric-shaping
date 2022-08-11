function C = RelabelNDQAM(C)

% check inputs
[M,D] = size(C);

% spread bits over dimensions
mtot = log2(M);
m = floor(mtot/D)*ones(1,D);

mlow = mtot-sum(m);
m(end-mlow+1:end) = m(end-mlow+1:end)+1;

% recursively map the constellation for every dimension
map = maprecursive(C,m);

I = eye(M);
if abs(det(I(:,map+1))) ~=1
    error('Did not find a map')
end

C(map+1,:) = C;
end

function map = maprecursive(x,m)
if size(x,2)==1
    % if the constellation is 1D, assign gray_mapped labels
    [~,argsort] = sort(x);
    map(argsort) = graymap(2^m(1));
else
    % if the constellation is ND
    map = zeros(2^sum(m),1);
    map_dim = graymap(2^m(1));
    % sort by the first column
    [~,argsort] = sortrows(x);
    % label the slices for the remainder of the dimensions individually
    for i=1:2^m(1)
        idx= argsort( (i-1)*2^sum(m(2:end)) + (1:2^sum(m(2:end))) );
        map(idx) = 2^sum(m(2:end))*map_dim(i) + maprecursive(x(idx,2:end),m(2:end));
    end
end
end

function mapping = graymap(order)
    j = int32((0:order-1)');
    mapping = cast(bitxor(j,bitshift(j,-1)),'like',order);
end