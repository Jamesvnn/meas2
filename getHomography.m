function H = getHomography(ms, M, eps_lim, max_iter)
    % ms - Sensor coordinates for every chessboard
    % M  - model plate
    H = initHomography(ms, M);

    N = size(M, 1);
    nImg = size(ms, 1); % number of captured images
    
    nrm = 1;
    mu = 0.001;
    for i = 1:nImg % i = 1:nImage
        mi = ms{i};
        Hi = H{i};
        
        H_iter = [Hi(1, 1), Hi(1, 2), Hi(1,3), Hi(2, 1), Hi(2, 2), Hi(2,3), Hi(3, 1), Hi(3, 2), Hi(3,3)]';
        
        k = 0;
        while 1
            k = k + 1;
            dJdH = [];
            difs  = [];
            for j = 1:N % calculation of epsilon & Jacobian
                tmp_tpX = Hi(1, 1)*M(j, 1) + Hi(1, 2)*M(j, 2) + Hi(1, 3);
                tmp_tpY = Hi(2, 1)*M(j, 1) + Hi(2, 2)*M(j, 2) + Hi(2, 3);
                tmp_bot = Hi(3, 1)*M(j, 1) + Hi(3, 2)*M(j, 2) + Hi(3, 3);
                
                difs = [difs; mi(j, 1) - tmp_tpX/tmp_bot; mi(j, 2) - tmp_tpY/tmp_bot];
            
                DJDH11 = -M(j, 1)/tmp_bot;
                DJDH12 = -M(j, 2)/tmp_bot;
                DJDH13 =       -1/tmp_bot;
                DJDH31 = tmp_tpX*M(j, 1)/tmp_bot^2;
                DJDH32 = tmp_tpX*M(j, 2)/tmp_bot^2;
                DJDH33 = tmp_tpX        /tmp_bot^2;
                
                dJdH = [dJdH; DJDH11 DJDH12 DJDH13 0 0 0 DJDH31 DJDH32 DJDH33];
                
                DJDH21 = -M(j, 1)/tmp_bot;
                DJDH22 = -M(j, 2)/tmp_bot;
                DJDH23 =       -1/tmp_bot;
                DJDH31 = tmp_tpY*M(j, 1)/tmp_bot^2;
                DJDH32 = tmp_tpY*M(j, 2)/tmp_bot^2;
                DJDH33 = tmp_tpY        /tmp_bot^2;
                
                dJdH = [dJdH; 0 0 0 DJDH21 DJDH22 DJDH23 DJDH31 DJDH32 DJDH33];
            end
            
            dHi = inv(dJdH'*dJdH + mu*ones(9))*dJdH'*difs;
            
%             disp(num2str(max(abs(dHi))));
            H_iter = H_iter - dHi;
            Hi = [
                    H_iter(1) H_iter(2) H_iter(3)
                    H_iter(4) H_iter(5) H_iter(6)
                    H_iter(7) H_iter(8) H_iter(9)                
                ];
            if (max(abs(dHi)) < eps_lim) | k > max_iter
                H{i} = Hi/Hi(3, 3);
                disp(['    Converged Homography - ', num2str(i), '(Error limit = ', num2str(max(abs(dHi))), ')'])
                break
            end
        end
    end   
end

function normal_data = normalMatrix(data)
    x = data(:, 1);
    y = data(:, 2);

    N = size(data, 1);

    x_mean = mean(x);
    y_mean = mean(y);
    x_var  = var(x)*(N-1)/N;
    y_var  = var(y)*(N-1)/N;
    
%     # Form rescaling matrix so that data points will lie
%     # sqrt(2) from the origin on average.
    s_x = sqrt(2. / x_var);
    s_y = sqrt(2. / y_var);
    
    normal_data = [
        s_x,  0., -s_x * x_mean;
        0., s_y, -s_y * y_mean;
        0.,  0.,            1.];

end

function H0 = initHomography(ms, M)
    % ms - Sensor coordinates for every chessboard
    % M  - model plate
    N = size(M,  1); % number of sensors
    nImg = size(ms, 1); % number of captured images

    H0 = {};
    norm_M = normalMatrix(M);
    new_M  = M*norm_M';
    for i = 1:nImg
        mi = ms{i};
        norm_mi = normalMatrix(mi);
        new_mi  = mi*norm_mi';
        
        L  = zeros(2*N, 9);
        L(1:N, 1:3) = new_M;
        L(1:N, 7:9) = -repmat(new_mi(:,1), 1, 3).*new_M();
        
        L(N+1:2*N, 4:6) = new_M;
        L(N+1:2*N, 7:9) = -repmat(new_mi(:,2), 1, 3).*new_M;
%         [V, D] = eig(L'*L);
%         
%         H0i = V(:,1);
        [U, S, V] = svd(L);
        [val, idx] = min(max(S));

        H0i = V(:, idx);
        H0i = reshape(H0i, 3, 3);
        % denormalize
        H0{i, 1} = inv(norm_mi)*H0i'*norm_M;
    end
end
