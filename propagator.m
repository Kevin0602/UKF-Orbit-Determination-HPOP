function  X_next= propagator( X_now,param )

X = Ephemeris(X_now(1:6,:), 1,param );
X_next = X(end,2:end)';
% if size(X,1) == 15 || size(X,1) == 14
%         X_next(1:6,:) = X_next(1:6,:) + X(7:14,:);
% end

end

