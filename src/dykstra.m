function x = dykstra(r,Pd,Pc,Niter)
% Approximates the solution to \min_x 0.5*||x - r||² s.t. x \in D \cap C
% (Dykstra algorithm).
%%
% Code : Pierre-Antoine Thouvenin, January 16th 2015.
%%
%-------------------------------------------------------------------------%
% Input: 
% > r       data to be projected;
% > Pd,Pc   function handles of the involved projectors;
%           (onto D and C respectively);
% > Niter   number of iterations.
%
% Output:
% < x       approximate projection onto D \cap C.
%-------------------------------------------------------------------------%
%%
p = zeros(size(r));
q = p;
x = r;

for k = 1:Niter
    y = Pd(x + p);
    p = x + p - y;
    x = Pc(y + q);
    q = y + q - x;
end

end
