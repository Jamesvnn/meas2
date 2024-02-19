function distorion = getDistortion(ms, M, intrinsic, extrinsic)
    % this code is based on Matlab's camera calibration expression
    Xc = []; % in camera coordinate
    for i = 1:length(extrinsic)
        ext_i = extrinsic{i}(:, [1, 2, 4]);
        for j = 1:size(M, 1)
            Xw_j = M(j, :);
            Xc_i = ext_i*Xw_j';
            
            Xc = [Xc; Xc_i'];
            
        end
    end
    
    x = [Xc(:, 1:2)]/Xc(3); % Undistorted pixel locations.
                            % x and y are in normalized image coordinates.
                            % Normalized image coordinates are calculated from pixel coordinates 
                            % by translating to the optical center and dividing by the focal length in pixels.
                            % Thus, x and y are dimensionless. 
    
    xd = [];
    A_inv = inv(intrinsic);
    for i = 1:length(ms)
        mi = ms{i};
        for j = 1:size(M, 1)
            mij = mi(j, :);
            xd = [xd; [A_inv*mij']'];
        end
    end
    
    xd = xd(:, 1:2); % distorted coordinates
    
    r2 = x(:, 1).^2 + x(:, 2).^2;
    r4 = r2.^2;
    
    D = [
        x(:, 1).*r2, x(:, 1).*r4, 2*x(:, 1).*x(:, 2), r2 + 2*x(:, 1).^2
        x(:, 2).*r2, x(:, 2).*r4, r2 + 2*x(:, 2).^2,  2*x(:, 1).*x(:, 2)
        ];
    d = (xd - x)';
    d = d(:);
    
%     D = D(:, 1:2);
    distorion = inv(D'*D)*D'*d;
end
