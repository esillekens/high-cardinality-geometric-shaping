function X = RelabelAPSK_spectral(X)


[M,D] = size(X);

% Only implemented for 2D
assert(D==2, 'Error, only implemented for 2D');


% bits
m = log2(M);

mA = floor(m/2);
mP = ceil(m/2);

% ensure more phase than amplitude
if mA==mP
    mA = mA-1;
    mP = mP+1;
end

% shift balance if single quadrant
if sum(all(X>=0))==2
    mA = mA+1;
    mP = mP-1;
end

% order
MP = 2^mP;
MA = 2^mA;

% phase and amplitude
phi = atan2(X(:,2),X(:,1));
amp = sqrt(sum(X.^2,2));

% prepare mapping
mapP = MA*graymap(MP);
mapA = graymap(MA);
mapOut = zeros(M,1);

% Use spectral clustering to group constellation points
idx_cluster = spectral_clustering(Polar2Cart(amp,mod(phi, pi/2/MP)), MA);

for i = 1:MA
    % select the i-th cluster
    idx = find(idx_cluster==i);
    % argsort the slice according to Phase
    [~,argsortA] = sort(phi(idx)); 
    % combine the bit labels
    mapOut(idx(argsortA)) = mapP + mapA(i);
end

% relabel the constellation
X(mapOut+1,:) = X;

end

function mapping = graymap(order)
    % generate graymapped labels
    j = int32((0:order-1)');
    mapping = cast(bitxor(j,bitshift(j,-1)),'like',order);
end

function idx = spectral_clustering(X, k)

% number of points to be clustered
M = size(X,1);

% caculate a likeness measure
SNR = 12;
gamma = 10.^(SNR/10);

EucD = pdist2(X,X,'squaredeuclidean');

% construct affinity matrix
A = exp(-gamma*EucD) - eye(M);
D = sum(A,1);

% convert to normalised Lagrangian matrix
L = diag(D) - A; %lagrangian
L = 1./sqrt(D).' .* L .* 1./sqrt(D); %normalisation

% extract eigen vectors
[x, ~] = eig(L);
x = x(:,1:k);
vectors = 1./sqrt(D).' .* x;

% Use the Fiedler vector to order the points
[~,iii] = sort(vectors(:,2));

% create output indices
idx = zeros(M,1);
idx(iii) = kron((1:k).',ones(M/k,1));

end

function X = Polar2Cart(amp,phi)

X = amp.*[cos(phi),sin(phi)];

end