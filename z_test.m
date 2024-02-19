function errors = z_test(ms, M, A, K, P, Rhos, Trns)
    clear
    clc
    
    load z_test.mat
    
%     load z_test_1.mat
    
%     A = A';
    k1 = K(1);
    k2 = K(2);
    p1 = P(1);
    p2 = P(2);
    for i = 1:length(Rhos)
        ms_real = ms{i};
        ms_proj = [];
        
        rot = Rodrigue2Rotation(Rhos{i});
%         rot = Rhos{i}';
        trn = Trns{i};
        
        W = [rot, trn'];
        for j = 1:size(M, 1)
            Xc = W(:, [1 2 4])*M(j, :)';
            
            x = Xc(1:2)/Xc(3);
            
            r2 = x(1)^2 + x(2)^2;
            xd(1) = x(1)*(1 + k1*r2 + k2*r2^2) + 2*p1*x(1)*x(2) + p2*(r2 + 2*x(1)^2);
            xd(2) = x(2)*(1 + k1*r2 + k2*r2^2) + 2*p2*x(1)*x(2) + p1*(r2 + 2*x(2)^2);
            
            ms_proj = [ms_proj; [A*[xd 1]']'];
        end
        
        errs = abs(ms_real - ms_proj);
        disp(num2str(max(errs(:))))
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