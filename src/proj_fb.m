function x = proj_fb(y,c,r)
% Projection of y onto the Frobenius ball B(c,r), of radius r > 0 and 
% center c
% 
%%
% Code : Pierre-Antoine Thouvenin, May 16th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > y       data matrix;
% > c       center of the Frobenius ball
% > r       radius of the Frobenius ball.
%
% Output:
% < x   projection of y ontot B(c,r)
%-------------------------------------------------------------------------%
%%
u = y-c;
x = c + u*min([1,r/norm(u,'fro')]); 
end



