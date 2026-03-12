function [Er_xyz, Etheta_xyz, Ephi_xyz] = ...
    sph2xyz(theta, phi, Er, Etheta, Ephi)

%
st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);



% e_r
erx = st .* cp;
ery = st .* sp;
erz = ct;

% e_theta
etx = ct .* cp;
ety = ct .* sp;
etz = -st;

% e_phi
epx = -sp;
epy =  cp;
epz = zeros(size(phi));



Er_xyz(:,:,1) = Er .* erx;
Er_xyz(:,:,2) = Er .* ery;
Er_xyz(:,:,3) = Er .* erz;

Etheta_xyz(:,:,1) = Etheta .* etx;
Etheta_xyz(:,:,2) = Etheta .* ety;
Etheta_xyz(:,:,3) = Etheta .* etz;

Ephi_xyz(:,:,1) = Ephi .* epx;
Ephi_xyz(:,:,2) = Ephi .* epy;
Ephi_xyz(:,:,3) = Ephi .* epz;

end
