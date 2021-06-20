function [xn] = vecnorm3d(x)

    xn = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);

end

