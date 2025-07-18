% create model setup from TIFF image

function [units,D,Nz] = ModelFromImageDiscrete(filename,W,Nx)

% read in RGB image from TIFF file
img = double(imread(filename)) / 255;

if size(img, 3) > 3
    img = img(:, :, 1:3);
end

% get size of image
[p,q,~] = size(img);
disp(p);
disp(q);

% RGB values in tiff image
yellow  = [1,230/255,128/255];
orange =  [1, 127/255, 39/255];
red =     [237/255, 28/255, 36/255];
white =   [1,1,1];

% mapping pixels to discrete values
img_units = zeros(p,q); %2D matrix same size as image
for i = 1:p
    for j = 1:q
        pixel = squeeze(img(i, j, :))';  % safely extract RGB vector

        if norm(pixel - yellow) < 0.01
            img_units(i,j) = 1; % yellow
        elseif norm(pixel - orange) < 0.01
            img_units(i,j) = 2; % orange
        elseif norm(pixel - red) < 0.01
            img_units(i,j) = 3; % red
        else
            img_units(i,j) = 4; % white (default)
        end
    end
end


% interpolate from original dimensions to target model size
D  = W*p/q;
Nz = floor(Nx*p/q);

ho  = W/q;
xco = ho/2:ho:W-ho/2;
zco = ho/2:ho:D-ho/2;
[Xco,Zco] = meshgrid(xco,zco);

h   = W/Nx;
xc  = h/2:h:W-h/2;
zc  = h/2:h:D-h/2;
[Xc,Zc] = meshgrid(xc,zc);

imgi = round(interp2(Xco,Zco,img_units,Xc,Zc,'nearest'));

figure(1); clf
subplot(2,1,1)
imagesc(xco,zco,img_units); axis equal tight; colorbar
subplot(2,1,2)
imagesc(xc,zc,imgi); axis equal tight; colorbar

units = uint8(imgi);

end