function [A, K, P, Rhos, Trns] = unpackParam(p)
    nT = length(p);
    A  = [
        p(1), p(2), p(3)
        0,    p(4), p(5)
        0,    0,    1
        ];
    
    K = [p(6); p(7)];
    P = [p(8); p(9)];
    
    nImg = (nT - 9)/6;
    
    if mod(nT - 9, 6) ~= 0
        Error('get_Jac_eps(): nImg must be integer...')
    end

    Rhos = {};
    Trns = {};
    for iImg = 1:nImg
        idx = 9 + (iImg - 1)*6;

        Rho(1) = p(idx + 1);
        Rho(2) = p(idx + 2);
        Rho(3) = p(idx + 3);

        Trn(1) = p(idx + 4);
        Trn(2) = p(idx + 5);
        Trn(3) = p(idx + 6);
        
        Rhos{iImg, 1} = Rho;
        Trns{iImg, 1} = Trn;
    end
end