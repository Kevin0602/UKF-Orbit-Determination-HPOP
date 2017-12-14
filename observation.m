function HX= observation( X ,param)

Z = X(1:3,:)'-param(:,1:3);
%normlize
Z = Z./(sqrt(Z(:,1).^2+Z(:,2).^2+Z(:,3).^2)*[1 1 1]);
% if size(X,1) == 15
%       Z = Z + X(15,:);
% end
HX=Z';

end

