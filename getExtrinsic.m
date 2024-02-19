function extrinsic = getExtrinsic(intrinsic, H)
    extrinsic = {};
    for i = 1:length(H)
        Hi = H{i};
        h1 = Hi(:, 1);
        h2 = Hi(:, 2);
        h3 = Hi(:, 3);
        
        lambda = 1/norm(inv(intrinsic)*h1, 2);
        
        r1 = lambda*inv(intrinsic)*h1;
        r2 = lambda*inv(intrinsic)*h2;
        r3 = cross(r1, r2);
        t  = lambda*inv(intrinsic)*h3;
        
        R = [r1 r2 r3];
        [U, ~, V] = svd(R);
        R = U*V';
        
        extrinsic{i, 1} = [R, t];
    end
end
