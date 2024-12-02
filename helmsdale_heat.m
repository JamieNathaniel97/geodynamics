% create model setup from TIFF image

function [units,D,Nz] = ModelFromImage(filename,n_units,W,Nx)

% read in RGB image from TIFF file
img = double(importdata(filename));

% get size of image
[p,q,k] = size(img);

% identify units by clustering analysis
[Ic,Fc] = kmeans(reshape(img,p*q,k),9,'MaxIter',1e3,'Replicates',5);

% sort clusters for ascending SiO2 content
[~,isort]   = sort(Fc(:,1)-sum(Fc,2),'descend');
Ics = zeros(size(Ic));
for ic = 1:n_units
    Ics(Ic==isort(ic)) = ic;
end


imgc = reshape(Ics,p,q);

% interpolate from original dimensions to target model size
D  = W*p/q;
Nz = floor(Nx*p/q);

%Coordinate grid on original image
ho  = W/q;
xco = ho/2:ho:W-ho/2;
zco = ho/2:ho:D-ho/2;
[Xco,Zco] = meshgrid(xco,zco);

%Coordinate grid on the output / for code interpretation: use this as the
%model grid 
h   = W/Nx;
xc  = h/2:h:W-h/2;
zc  = h/2:h:D-h/2;
[Xc,Zc] = meshgrid(xc,zc);


%Set up Insulating Top and Bottom 
iz3 = [           1,1:Nz,Nz      ];  
iz5 = [        1, 1,1:Nz,Nz,Nz   ];
iz7 = [   1,   1, 1,1:Nz,Nz,Nz,Nz];



switch ADVN
    case 'WENO5'
        fxppos = weno5poly(fxmm ,  fxm, fxc, fxp, fxpp);     fxpneg = weno5poly(fxppp, fxpp, fxp, fxc, fxm );
        fxmpos = weno5poly(fxmmm, fxmm, fxm, fxc, fxp );     fxmneg = weno5poly( fxpp,  fxp, fxc, fxm, fxmm);

        fzppos = weno5poly(fzmm ,  fzm, fzc, fzp, fzpp);     fzpneg = weno5poly(fzppp, fzpp,fzp, fzc, fzm );
        fzmpos = weno5poly(fzmmm, fzmm, fzm, fzc, fzp );     fzmneg = weno5poly( fzpp, fzp, fzc, fzm, fzmm);
end





% Function to calculate diffusion rate
function [dTdt] = diffusion(f,k,h,ix,iz)
% calculate heat flux by diffusion

qx = - k .* diff(f(ix, :), 1, 2)/h;
qz = - k .* diff(f(iz, :), 1, 1)/h;

% calculate flux balance for rate of change
dTdt = - (diff(qx)/h+diff(qz)/h);


end


%*****  calculate numerical error norm
Err = norm(T - Ta,2)./norm(Ta,2);
disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error = ',num2str(Err)]);
disp(' ');

%plot the interpolated image
%imgi = interp2(Xco,Zco,imgc,Xc,Zc);
%figure(1); clf
%subplot(2,1,1)
%imagesc(xco,zco,imgc); axis equal tight; colorbar
%subplot(2,1,2)
%imagesc(xc,zc,imgi); axis equal tight; colorbar
%
%units = uint8(imgi);

end