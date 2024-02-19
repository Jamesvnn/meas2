clear
close all
clc

load mM.mat

disp('- Calculating Homographies...')
H  = getHomography(ms, M, 1e-20, 10);

disp('- Calculating Intrinsic parameters...')
[intrinsic, ~] = getIntrinsic(H);

disp('- Calculating Extrinsic parameters...')
extrinsic = getExtrinsic(intrinsic, H);

disp('- Calculating Distortion parameters...')
distortion = getDistortion(ms, M, intrinsic, extrinsic);

disp('- Adjusting all parameters finally...')
params = refinement(ms, M, intrinsic, extrinsic, distortion);
[A, K, P, Rhos, Trns] = unpackParam(params);

% errors = z_test(ms, M, A, K, P, Rhos, Trns);
% 
% K = distortion(1:2);
% P = distortion(3:4);
% A = intrinsic;

save z_test.mat ms M A K P Rhos Trns