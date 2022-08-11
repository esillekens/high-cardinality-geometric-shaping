function [GMI,dMIdx] = GMI_withgrad_nD(SNR,x)

[M,D] = size(x);

withgrad = nargout>=2;

% disp(withgrad)

%% correct for SNR
sigma_z = 10^(-SNR/20);

input_scaling = sqrt(mean(sum(x.^2,2)/(D/2)));

x = x./input_scaling; % unit power per 2D
x = x./sigma_z;

%% check for symmetry and remove mirror images
m = log2(M);
map = logical(de2bi(0:M-1,m));

pos_idx = ones(M,1,'logical');
sym_bits = zeros(m,1);
for d = 1:D
    [~,sym_bit] = max(abs(x(:,d).'*map));
    flip = ones(1,D);
    flip(d) = -1;
    ids = map(:,sym_bit)==0;
    if all(sum(abs(x(ids,:)-flip.*x(~ids,:)).^2,2)<eps(2*M/sigma_z))
        pos_idx(~ids) = 0;
        sym_bits(sym_bit) = d;
    end
end
pos_idx = find(pos_idx);
%%

% get GaussHermite Weights
L =10;
[xi,al] = GaussHermite(L);

% span over N dimensions
[z,alpha] = grid_nD(D,xi,al);

% Optionally only use large alpha 
% sel = alpha>=1e-6;
% z = z(sel,:);
% alpha = alpha(sel);

sumI = zeros(M,1);
sumP = zeros(M,m);
ddx = zeros(M,D);
temp_j = zeros(D,M-1);
temp_k = zeros(D,M/2-1);

% move to GPU if large enough
if (gpuDeviceCount>1)&&numel(x)>=2^7
        z = gpuArray(z);
        alpha = gpuArray(alpha);
        x = gpuArray(x);
        ddx = gpuArray(ddx);
        sumI = gpuArray(sumI);
        sumP = gpuArray(sumP);
        temp_j = gpuArray(temp_j);
        temp_k = gpuArray(temp_k);
end

for idx = 1:length(pos_idx)
    i = pos_idx(idx);
    dij = (x(i,:)-x((1:M)~=i,:)).';
    mapi = (map==map(i,:));
    mapi(i,:) = false;
    exp_n = exp(z*(-2*dij)-sum(dij.^2,1));

    sumJ = sum(exp_n,2);
    sumI(i) = alpha.' * log1p(sumJ);
    
    if withgrad
        den = exp_n./(1+sumJ);
        for d = 1:D
            temp_j(d,:) = -2*alpha.' *(den.*(z(:,d)+dij(d,:)));
        end

        ddx((1:M)~=i,:) = ddx((1:M)~=i,:) - m*temp_j.';
        ddx(i,:) = ddx(i,:) + m*sum(temp_j.');
    end
    
    for k=1:m
        exp_k = exp_n(:,mapi((1:M)~=i,k));
        sumJ = sum(exp_k,2);

        sumP(i,k) = alpha.'* log1p(sumJ);
        if withgrad
            den = exp_k./(1+sumJ);
            for d = 1:D
                temp_k(d,:) = -2*alpha.' *...
                    (den.*...
                        (z(:,d)+dij(d,mapi((1:M)~=i,k)))...
                    );
            end

            ddx(mapi(:,k),:) = ddx(mapi(:,k),:) + temp_k.';
            ddx(i,:) = ddx(i,:)-sum(temp_k.');
        end
    end
end
if withgrad
    for k=find(sym_bits~=0).'
        d = sym_bits(k);
        ids = map(:,k)==0;
        flip = ones(1,D);
        flip(d) = -1;
        ddx(ids,:) = ddx(ids,:) + flip.*ddx(~ids,:);
        ddx(~ids,:) = flip.*ddx(ids,:);
    end
end


Msum = numel(pos_idx);
GMI = log2(M) - m/Msum/pi^(D/2)*sum(sumI)/log(2)+1/Msum/pi^(D/2)*sum(sumP,'all')/log(2);
GMI = cast(gather(GMI),'like',SNR);

if withgrad
    dMIdx = -1/Msum/pi^(D/2)/sigma_z/input_scaling/log(2)*(ddx);
    dMIdx = cast(gather(dMIdx),'like',SNR);
end


end

function [z,alpha] = grid_nD(D,xi,al)

L = length(xi);

z_cell = cell(D,1);
alpha_cell = cell(D,1);

z_cell{1} = xi;
alpha_cell{1} = al;
for d = 2:D
    z_cell{d} = [repmat(z_cell{d-1},L,1),kron(xi,ones(L^(d-1),1))];
    alpha_cell{d} = repmat(alpha_cell{d-1},L,1).*kron(al,ones(L^(d-1),1));
end

z = z_cell{D};
alpha= alpha_cell{D};

end