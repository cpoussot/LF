function Hm = mirror(H)

% >> Mirror
[Z,P,G] = zpkdata(H);
nProj   = 0;
for ii = 1:length(P{1})
    % Mirror projection
    if real(P{1}(ii)) > 0
        P{1}(ii) = complex(-real(P{1}(ii)),imag(P{1}(ii)));
        nProj    = nProj + 1;
    end
    % Pole close to axis set to zero 
    if norm(P{1}(ii)) < 1e-3
        P{1}(ii) = 0;
    end
end
Hm      = ss(zpk(Z,P,G));
Hm.e    = eye(length(Hm.a));

% Hs{2} =  ss(zpk(Z,P,G));
% N(2)  = norm(squeeze(freqresp(Hs{2},W_shift) - Kstar_shift));
% Hs{3} = -ss(zpk(Z,P,G));
% N(3)  = norm(squeeze(freqresp(Hs{3},W_shift) - Kstar_shift));
% % Chose 
% [~,idx] = sort(N,'ascend');
% Hm      = Hs{idx(1)};
