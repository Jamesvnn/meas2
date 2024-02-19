function param = refinement(ms, M, intrinsic, extrinsic, distortion)
    eps = 1e-10;
    
    param = pack_param(intrinsic, extrinsic, distortion);
    
    iter = 0;
    mu = 1e-8;
    
    while 1
        iter = iter + 1;
        J     = [];
        e     = [];
        for i = 1:length(ms)
            mi = ms{i};
            [A, k1, k2, p1, p2, Rho, Trn] = getParam(param, i);
            
            for j = 1:size(M, 1)
                mij = mi(j, :);
                Xw  = M (j, :);
                
                uij = eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2);

                A_orig = A;
                A(1, 1) = A(1, 1) + eps;
                DJdA11 = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;

                A = A_orig;
                A(1, 2) = A(1, 2) + eps;
                DJdA12 = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;

                A = A_orig;
                A(1, 3) = A(1, 3) + eps;
                DJdA13 = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;

                A = A_orig;
                A(2, 2) = A(2, 2) + eps;
                DJdA22 = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;

                A = A_orig;
                A(2, 3) = A(2, 3) + eps;
                DJdA23 = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;

                A = A_orig;
                
                Rho_orig = Rho;
                Rho(1)   = Rho(1) + eps;
                DJdR1    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Rho      = Rho_orig;
                
                Rho(2)   = Rho(2) + eps;
                DJdR2    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Rho      = Rho_orig;
                
                Rho(3)   = Rho(3) + eps;
                DJdR3    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Rho      = Rho_orig;
                
                Trn_orig = Trn;
                Trn(1)   = Trn(1) + eps;
                DJdT1    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Trn      = Trn_orig;
                
                Trn(2)   = Trn(2) + eps;
                DJdT2    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Trn      = Trn_orig;
                
                Trn(3)   = Trn(3) + eps;
                DJdT3    = (eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2) - uij)/eps;
                Trn      = Trn_orig;

                DJdk1 = (eval_fullsystem(Xw, A, Rho, Trn, k1 + eps, k2      , p1      , p2      ) - uij)/eps;
                DJdk2 = (eval_fullsystem(Xw, A, Rho, Trn, k1      , k2 + eps, p1      , p2      ) - uij)/eps;
                DJdp1 = (eval_fullsystem(Xw, A, Rho, Trn, k1      , k2      , p1 + eps, p2      ) - uij)/eps;
                DJdp2 = (eval_fullsystem(Xw, A, Rho, Trn, k1      , k2      , p1      , p2 + eps) - uij)/eps;

                Jij = [
                    DJdA11(1:2)', DJdA12(1:2)', DJdA13(1:2)', DJdA22(1:2)', DJdA23(1:2)',...
                    DJdk1(1:2)' , DJdk2(1:2)' , DJdp1(1:2)' , DJdp2(1:2)', ...
                    zeros(2, (i-1)*6),...
                    DJdR1(1:2)', DJdR2(1:2)', DJdR3(1:2)', DJdT1(1:2)', DJdT2(1:2)', DJdT3(1:2)',...
                    zeros(2, (length(ms)-i)*6)
                    ];
                

                J = [J; Jij];
                
                eij = uij - mij;
                e   = [e; eij(1:2)'];
            end
        end
        
        d_param = inv(J'*J + mu*eye(size(J, 2)))*J'*e;
        param = param - d_param;
        
        e_norm = max(abs(e(:)));
        
        disp(['Iter = ', num2str(iter), ', Norm = ', num2str(max(abs(d_param)))])
        
        if iter > 5
            mu = 10*mu;
        end
        if max(abs(d_param)) < eps | iter > 100
            disp(['Iter = ', num2str(iter), ', Norm = ', num2str(max(abs(d_param)))])
            break
        end
        
        e_norm_old = e_norm;
    end
end

function u = eval_fullsystem(Xw, A, Rho, Trn, k1, k2, p1, p2)
    W  = Rodrigue2Rotation(Rho);
    
    W  = [W Trn'];
    Xc = W(:, [1, 2, 4])*Xw';
    x  = Xc(1:2)/Xc(3);
    
    r2 = x(1)^2 + x(2)^2;
    xd(1) = x(1)*(1 + k1*r2 + k2*r2^2) + 2*p1*x(1)*x(2)     + p2*(r2 + 2*x(1)^2);
    xd(2) = x(2)*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*x(2)^2) + 2*p2*x(1)*x(2);
    
    u = A*[xd 1]';
    
    u = u';
end

function rho = Rotation2Rodrigue(R)
    R11 = R(1, 1);
    R12 = R(1, 2);
    R13 = R(1, 3);
    R21 = R(2, 1);
    R22 = R(2, 2);
    R23 = R(2, 3);
    R31 = R(3, 1);
    R32 = R(3, 2);
    R33 = R(3, 3);
    
    p = 0.5*[R32 - R23, R13 - R31, R21 - R12];
    C = (trace(R) - 1)/2;
    
    if norm(p, 2) == 0 & C == 1
        rho(1:3) = 0;
    elseif norm(p, 2) == 0 & C == -1
        Rp = R + eye(3);
        
        nrm(1) = norm(Rp(:, 1), 2);
        nrm(2) = norm(Rp(:, 2), 2);
        nrm(3) = norm(Rp(:, 3), 2);

        [~, idx] = max(nrm);
        
        v = Rp(:, 3);
        
        u = v/norm(v, 2);
        rho = pi*s;
    elseif norm(p, 2) ~= 0
        u = p/norm(p, 2);
        
        if C > 0
            theta = atan(norm(p, 2)/C);
        elseif C < 0
            theta = pi + atan(norm(p, 2)/C);
        elseif C == 0 & y > 0
            theta = pi/2;
        elseif C == 0 & y < 0
            theta = -pi/2;
        else
            Error('Rotation2Rodrigue');
        end
        rho = theta*u;
    else
        rho = [];
    end
end

function s = S(x)
    x0 = x(1);
    x1 = x(2);
    x2 = x(3);
    
    if(x0 < 0) & (x0 == 0 | x1 < 0) & ((x0 == 0 & x1 == 0) | x2 < 0)
        s = -x;
    else
        s = x;
    end
end

function rot = Rodrigue2Rotation(rho)
    theta = norm(rho, 2);
    rho = rho/norm(rho, 2);
    
    W = [
         0,      -rho(3),   rho(2)
         rho(3),  0,       -rho(1)
        -rho(2),  rho(1),   0
        ];
    rot = eye(3) + W*sin(theta) + W*W*(1 - cos(theta));
end

function p = pack_param(intrinsic, extrinsic, distortion)
    A  = intrinsic;
    k1 = distortion(1);
    k2 = distortion(2);
    p1 = distortion(3);
    p2 = distortion(4);

    p = [A(1, 1), A(1, 2), A(1, 3), A(2, 2), A(2, 3), k1, k2, p1, p2];
    for i = 1:length(extrinsic)
        W  = extrinsic{i};
        
        Rho = Rotation2Rodrigue(W(:, 1:3));
        
        if (length(Rho) == 0)
                Error('Wrong Rodrigue Transformation...')
        end
        
        Trn = W(:, 4);
        p = [p, Rho(1), Rho(2), Rho(3), Trn(1), Trn(2), Trn(3)];
    end
    
    p = p';
end

function [A, k1, k2, p1, p2, Rho, Trn] = getParam(p, iImg)
    nT = length(p);
    A  = [
        p(1), p(2), p(3)
        0,    p(4), p(5)
        0,    0,    1
        ];
    
    k1 = p(6);
    k2 = p(7);
    p1 = p(8);
    p2 = p(9);
    
    nImg = (nT - 9)/6;
    
    if iImg > nImg
        Error('get_Jac_eps():  iImg > nImg...')
    end
    
    if mod(nT - 9, 6) ~= 0
        Error('get_Jac_eps(): nImg must be integer...')
    end

    idx = 9 + (iImg - 1)*6;
    
    Rho(1) = p(idx + 1);
    Rho(2) = p(idx + 2);
    Rho(3) = p(idx + 3);
    
    Trn(1) = p(idx + 4);
    Trn(2) = p(idx + 5);
    Trn(3) = p(idx + 6);
end