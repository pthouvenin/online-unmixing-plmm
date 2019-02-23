function indices = outlier_patch(h_patch,w_patch,H,W,order)
% Calcul des coordonnées des points du patch (ordre lexicographique
% sous-jacent) à partir des coordonnées pixeliques absolues.
if max(w_patch(:)) > W || max(h_patch(:)) > H
    error('The pixel coordinates exceed the image size.');
end

if order
    % ordre lexicographique
    h_patch = h_patch';
    indices = (h_patch(:,ones(1,length(w_patch))) - 1)*W + w_patch(ones(1,length(h_patch)),:);
    indices = indices(:);
else
    % ordre non-lexicographique
    w_patch = w_patch';
    indices = (w_patch(:,ones(1,length(h_patch))) - 1)*H + h_patch(ones(1,length(w_patch)),:);
    indices = indices(:);
end