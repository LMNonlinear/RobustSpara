function Jinv=pinv_ridge(J,lambda)
% pinv_ridge=@(J,lambda) (J'*J+lambda*eye(size(J,2)))\J';
% Jinv=(J'*J+lambda*eye(size(J,2)))\J';
% 
N=size(J,1);
% lambda=lambda*trace(J'*J)/N;
lambda=lambda*sum(J.*J,'all')/N;
Jinv=(J'*J+lambda*eye(size(J,2)))\J';
end