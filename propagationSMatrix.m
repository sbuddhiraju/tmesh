function S_matrix = propagationSMatrix( thickness_list, M_matrix, f_matrix, numOfG, startLayer)
%This function will return the S_matrix at (kx,ky)
% I_matrix, f_matrix:                        matrices from fluxAtKxKy
% numOfG:                                    number of G chosen
% startLayer:                                the layer that propagation starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:    Kaifeng (Francis) Chen
% Copyright: Shanhui Fan's Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PhysicsConst;

N = 2*numOfG + 1;
numOfLayer = length(thickness_list);
r1 = 1:2*N;
r2 = 2*N+1:4*N;
%% definining scattering matrices
S_matrix = cell(numOfLayer,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numOfLayer
    S_matrix{i} = eye(4*N);
end

%% calculating the S matrix
%% propagating below sourceLayer
for i = startLayer:-1:2
    I_matrix = M_matrix{i} \ M_matrix{i-1};
    %% important step, reverse I according to the order of a and b
    I_matrix = [I_matrix(r2,r2), I_matrix(r2,r1) ; I_matrix(r1,r2), I_matrix(r1,r1)];
    S_matrix{i-1}(r1,r1) = (I_matrix(r1,r1) - f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r1))\f_matrix{i}*S_matrix{i}(r1,r1);
    S_matrix{i-1}(r1,r2) = (I_matrix(r1,r1) - f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r1))\...,
        (f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r2) - I_matrix(r1,r2))*f_matrix{i-1};
    S_matrix{i-1}(r2,r1) = S_matrix{i}(r2,r1) + S_matrix{i}(r2,r2)*I_matrix(r2,r1)*S_matrix{i-1}(r1,r1);
    S_matrix{i-1}(r2,r2) = S_matrix{i}(r2,r2)*(I_matrix(r2,r1)*S_matrix{i-1}(r1,r2)+I_matrix(r2,r2)*f_matrix{i-1});
end

%% propagating after sourceLayer
for i = startLayer:numOfLayer-1
    I_matrix = M_matrix{i} \ M_matrix{i+1};
    S_matrix{i+1}(r1,r1) = (I_matrix(r1,r1) - f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r1))\f_matrix{i}*S_matrix{i}(r1,r1);
    S_matrix{i+1}(r1,r2) = (I_matrix(r1,r1) - f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r1))\...,
        (f_matrix{i}*S_matrix{i}(r1,r2)*I_matrix(r2,r2) - I_matrix(r1,r2))*f_matrix{i+1};
    S_matrix{i+1}(r2,r1) = S_matrix{i}(r2,r1) + S_matrix{i}(r2,r2)*I_matrix(r2,r1)*S_matrix{i+1}(r1,r1);
    S_matrix{i+1}(r2,r2) = S_matrix{i}(r2,r2)*(I_matrix(r2,r1)*S_matrix{i+1}(r1,r2)+I_matrix(r2,r2)*f_matrix{i+1});
end

end


