clear
close all
clc

images = imageDatastore(fullfile(toolboxdir('vision'), 'visiondata', ...
                                    'calibration', 'mono'));

[imagePoints, boardSize] = detectCheckerboardPoints(images.Files);

squareSize  = 29;
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

I                   = readimage(images,1); 
imageSize           = [size(I, 1), size(I, 2)];

tic
[params, ~, errors] = estimateCameraParameters(imagePoints, worldPoints, ...
                                  'ImageSize', imageSize,...
                                  'EstimateTangentialDistortion', true,...
                                  'EstimateSkew' , true);
toc
