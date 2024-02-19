function [intrinsic, lambda] = getIntrinsic(H)
    nImg = length(H); % number of captured images
    V = [];
    for i = 1:nImg
        Hi = H{i}';

        V11 = getV_ij(Hi, 1, 1);
        V12 = getV_ij(Hi, 1, 2);
        V22 = getV_ij(Hi, 2, 2);

        V = [V; V12; V11 - V22];
    end

    [~, S, R]  = svd(V);
    [~, idx] = min(max(S));

    b = R(:, idx);
    
    B11 = b(1);
    B12 = b(2);
    B22 = b(3);
    B13 = b(4);
    B23 = b(5);
    B33 = b(6);
    
    v0     = (B12*B13 - B11*B23)/(B11*B22 - B12^2);
    lambda = B33 - (B13^2 + v0*(B12*B13 - B11*B23))/B11;
    alpha  = sqrt(lambda/B11);
    beta   = sqrt(lambda*B11/(B11*B22 - B12^2));
    gamma  = -B12*alpha^2*beta/lambda;
    u0     = gamma*v0/beta - B13*alpha^2/lambda;
    
    intrinsic = [
        alpha gamma u0
        0     beta  v0
        0     0     1
        ];
end

function Vij = getV_ij(H, i, j)
        Vij    = zeros(1, 6);
        Vij(1) = H(i, 1)*H(j, 1);
        Vij(2) = H(i, 1)*H(j, 2) + H(i, 2)*H(j, 1);
        Vij(3) = H(i, 2)*H(j, 2);
        Vij(4) = H(i, 3)*H(j, 1) + H(i, 1)*H(j, 3);
        Vij(5) = H(i, 3)*H(j, 2) + H(i, 2)*H(j, 3);
        Vij(6) = H(i, 3)*H(j, 3);
    end