function C = MapASPSK(C)

[M,D] = size(C);

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
if sum(all(C>=0))==2
    mA = mA+1;
    mP = mP-1;
end

% order
MP = 2^mP;
MA = 2^mA;

% phase and amplitude
phi = atan2(C(:,2),C(:,1));
amp = sqrt(sum(C.^2,2));

% prepare mapping
mapP = MA*graymap(MP);
mapA = graymap(MA);
mapOut = zeros(M,1);

% sort by phase
[~,argsortP] = sort(phi);

for i = 1:MP
    idx = argsortP((i-1)*MA + (1:MA)); % idx in the pizza slice
    [~,argsortA] = sort(amp(idx)); % argsort the slice according to amplitude
    mapOut(idx(argsortA)) = mapA + mapP(i);  
end

C(mapOut+1,:) = C;
end

function mapping = graymap(order)
    j = int32((0:order-1)');
    mapping = cast(bitxor(j,bitshift(j,-1)),'like',order);
end