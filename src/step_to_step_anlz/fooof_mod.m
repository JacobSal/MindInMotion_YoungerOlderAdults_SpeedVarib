function psd = fooof_mod(F,a,b,k,A,W,C)
    % a = aperiodic exponent
    % b = aperiodic offset
    % k = knee parameter
    % A = vector (nx1) of the power at peak n
    % W = vector (nx1) of the standard deviation at peak n
    % C = vector (nx1) of the center freqeuncy at peak n
    L = b - log10(k + F.^a);
    G = zeros(length(A),length(F));
    for n = 1:length(A)
        G(n,:) = A(n) .* exp((-(F-C(n)).^2)./(2*W(n).^2));
    end
    % G(n,:) = A * exp((-(F-C).^2)./(2*W.^2));
    psd = L + sum(G,1);
end