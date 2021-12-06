function [XX] = generate_rnd_pts_in_ss(Ndim,r1,Npoints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X = randn(Npoints,Ndim); % A normal dist is symmetrical
X = r1*X./sqrt(sum(X.^2,2)); % project to the surface of the Ndim-sphere
R = nthroot(rand(Npoints,1),Ndim);
% Combine
XX = X.*R;
end

