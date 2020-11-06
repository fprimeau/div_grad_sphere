function [S] = div_grad(lat,lon);
% [div, grad, del2] = div_grad(lat,lon);
% output:
%    div:    2-d divergence operator on the sphere
%    grad:   2-d divergence operator on the sphere
% input:
%    lat: vector of equally spaced latitudes in degrees
%    lon: vector of equally spaced longitudes in degrees
%    lat should not include either the north pole or the south pole.
% FWP 11/5/2020 created but not tested
% If someone does test it and confirms that it works or doesn't please let me know. THNX!
% FWP fprimeau@uci.edu     
  
    %% Sample usage:                                                                  
    % Find the first few eigenmodes of Laplace's equation on the surface of a the Earth
    %%                                                                                
    % dlat = 1; dlon = 1; 
    % lon = dlon/2:dlon:360-(dlon/2);
    % lat = (-90+dlat/2):dlat:(90-dlat/2);
    % S = div_grad(lat,lon);
    % [V,D] = eigs(S.del2,64,'smallestabs');
    % Q = S.DX+nan; 
    % LAM = S.LAM; PHI =  S.PHI; a = S.a;
    % x = a*cos(PHI).*cos(LAM);
    % y = a*cos(PHI).*sin(LAM);
    % z = a*sin(PHI);
    % figure(1); clf
    % for i = 1:9
    %   subplot(3,3,i)
    %   Q(:) = V(:,i);
    %   surf(x,y,z,Q); axis equal; shading flat
    % end
    
    %
    ny = length(lat);  nx = length(lon);  n = nx*ny;
    dlat = lat(2)-lat(1);     dphi = dlat*pi/180;
    dlon = lon(2)-lon(1);     dlam = dlon*pi/180;   
    % compute the geometric average of the polar and equatorial radii
    b = 6378e3; % equatorial radius [m]
    c = 6356e3; % polar radius [m]
    a = (b*b*c)^(1/3); % geometrical mean radius [m]
    % make a 2d mesh of points
    [LON,LAT] = meshgrid(lon,lat);          % [deg]
    PHI = LAT*pi/180;  LAM = LON*pi/180;    % [rad]
    DY = a*dphi+0*LAT; DX = a*cos(PHI)*dlam;% [m]

    % make doubly periodic shift operators
    ii = zeros(ny,nx);
    ii(:) = 1:n;
    I = speye(n);
    ie = ii(:,[2:end,1]); 
    in = ii([2:end,end],:); 
    IE = I(ie,:); IN = I(in,:); 
    % make a gradient operator
    d0 = @(x)  spdiags(x(:),0,length(x(:)),length(x(:))); % sparse diag function
    grad = [d0(DX(:))\(IE-I);...
            d0(DY(:))\(IN-I)];
    div = -grad.';
    del2 = div*grad;
    % put the output in a convenient stucture
    S.a = a;
    S.DX = DX;
    S.DY = DY;
    S.DA = DX.*DY;
    S.PHI = PHI;
    S.LAM = LAM;
    S.div = div;
    S.grad = grad;
    S.del2 = del2;
end

