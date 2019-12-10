function [flux, work] = fluxAtKxKy( omega, thickness_list, kx, ky_list,dielectricMatrix, dielectricMatrixInverse, dielectric_im_Matrix, numOfG, sourceLayer,targetLayer,modLayer, period,angle)
%This function will return the flux from the source layer to target layer
%This function will also perform the z integral
% kx                                         the parellel kx vector
% ky_list                                    parellel ky vector list, normalized
% omega:                                     angular frequency
% dielectricMatrix, dielectricMatrixInverse: the matrices used to calculate
%                                            the eigenmodes
% dielectric_im_Matrix:                      Imaginary part of the epsilon
% numOfG:                                    number of G chosen
% numOfLayer:                                number of layers in the system
% sourceLayer:                               place source in this layer
% modLayer:                                  where work is computed
% targetLayer:                               where the flux is measured
% period:                                    temporal periodicity
% angle:                                     'all' or 'normal'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tMESH: Siddharth Buddhiraju
% Based on MESH by Kaifeng Chen
% Copyright: Shanhui Fan's Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PhysicsConst;

N = 2*numOfG + 1; %number of modulation bands
Omega = 2 * pi / period;
W = omega * eye(N) + diag(-numOfG:numOfG)*Omega/c;

numOfLayer = length(thickness_list);
ky_list = ky_list * omega;
kx_matrix = eye(N)*kx;
flux = zeros(N,length(ky_list));
work = zeros(1,length(ky_list));

%% Setup matrices based on MESH paper
%For variable definitions in this function, please refer to the MESH paper

E_matrix = cell(numOfLayer,1);
H_matrix = cell(numOfLayer,1);
T_matrix = cell(numOfLayer,1);
M_matrix = cell(numOfLayer,1);

r1 = 1:2*N;
r2 = 2*N+1:4*N;

%% Obtain E and H for each layer
for i=1:numOfLayer
    E_matrix{i} = blkdiag(dielectricMatrix{i},dielectricMatrix{i});
    H_matrix{i} = blkdiag(dielectricMatrixInverse{i},dielectricMatrixInverse{i});
end

%% Define eigenvalue array and eigenmode matrix
q = zeros(N,1);
eigenModes = cell(numOfLayer,1);
qValues = cell(numOfLayer,1);
f_matrix = cell(numOfLayer, 1);

%% Loop over ky_list
for i=1:length(ky_list)
    
    ky = ky_list(i);
    K_matrix = [kx_matrix*kx_matrix, ky * kx_matrix;...
        ky * kx_matrix,ky^2 * eye(N)]; %see MESH paper
    
    for j=1:numOfLayer
        T_matrix{j} = [ky^2*dielectricMatrixInverse{j},-ky*dielectricMatrixInverse{j}*kx_matrix;...,
            -ky*kx_matrix*dielectricMatrixInverse{j},kx_matrix*dielectricMatrixInverse{j}*kx_matrix];
        
        eigenMatrix = E_matrix{j}*(kron(eye(2),W^2) - T_matrix{j})-K_matrix;
        [ eigenModes{j},qSquare] = eig(eigenMatrix); %solve eigenmatrix
        
        for k = 1:2*N
            q(k) = sqrt(qSquare(k,k)); %obtain z- propagation constants
            if imag(q(k)) > 0
                q(k) = -q(k);
            end
            
        end
        
        qValues{j} = diag(q);
        
        if thickness_list(j) == inf
            f_matrix{j} = eye(2*N);
        else
            f_matrix{j} = diag(exp(-1i*q*thickness_list(j)));
        end
        
        % Create M matrices, see MESH paper for definition
        M_matrix{j} = [(kron(eye(2),W^2) - T_matrix{j})*eigenModes{j}/qValues{j},...
            -(kron(eye(2),W^2) - T_matrix{j})*eigenModes{j}/qValues{j};...,
            kron(eye(2),W)*eigenModes{j}, kron(eye(2),W)*eigenModes{j}];
        
    end
    
    %% Begin calculation of flux
    S_matrix_target = propagationSMatrix( thickness_list, M_matrix, f_matrix, numOfG, targetLayer);
    
    fdt = zeros(N); %matrix to place thermal source in central band
    midpoint = (N+1)/2; %the central frequency order
    fdt(midpoint,midpoint) = 1; %source only present in central band
    
    grandImaginaryMatrix = blkdiag(fdt*dielectric_im_Matrix{sourceLayer}*fdt,...
        fdt*dielectric_im_Matrix{sourceLayer}*fdt, fdt*dielectric_im_Matrix{sourceLayer}*fdt);
    [q_R, q_L] = meshgrid(diag(qValues{sourceLayer}));
    
    source= [zeros(N), zeros(N), -ky*(dielectricMatrixInverse{sourceLayer}/W);...,
        zeros(N), zeros(N), kx_matrix*(dielectricMatrixInverse{sourceLayer}/W);...,
        zeros(N), eye(N), zeros(N);...,
        -eye(N), zeros(N), zeros(N);]; %see MESH paper for source
    
    %% Obtain EM field in target layer
    new_f_matrix = f_matrix;
    new_f_matrix{sourceLayer} = eye(2*N);
    S_matrix = propagationSMatrix( thickness_list, M_matrix, new_f_matrix, numOfG, sourceLayer);
    targetFields = M_matrix{sourceLayer} \ source; %inv(Ms)*P, Eq. (43, 51)
    
    %% Calculate P_1 and P_2
    P_1 = (eye(2*N) - S_matrix{targetLayer}(r1,r2) * S_matrix_target{numOfLayer}(r2,r1)) \ S_matrix{targetLayer}(r1,r1);
    P_2 = S_matrix_target{numOfLayer}(r2,r1) * P_1;
    %these are D_1 and D_2 in MESH paper
    
    %% Calculate Q_1 and Q_2
    Q_1 = eye(2*N) - f_matrix{sourceLayer} * S_matrix{1}(r2,r1) * f_matrix{sourceLayer} * S_matrix{numOfLayer}(r2,r1);
    Q_2 = -f_matrix{sourceLayer} * S_matrix{1}(r2,r1);
    % these are Kap_1 and Kap_2 in MESH paper
    
    %% Calculate R
    R = M_matrix{targetLayer} * [f_matrix{targetLayer} * P_1 ; P_2] / Q_1 * [eye(2*N), Q_2];
    % this is R in MESH paper (Eq. 43)
    
    %% Calculate Y*
    if thickness_list(sourceLayer) == inf
        integralMatrix_self = 1./(1i*(q_L-conj(q_R)));
        integralMatrix_mutual = zeros(size(q_L));
    else
        integralMatrix_self = (1-exp(-1i*(q_L-conj(q_R))*thickness_list(sourceLayer)))./(1i*(q_L-conj(q_R)));
        integralMatrix_mutual = (exp(1i*conj(q_R)*thickness_list(sourceLayer))-exp(-1i*q_L*thickness_list(sourceLayer)))./(1i*(q_L+conj(q_R)));
    end
    
    %% Calculate JZJ.Y*
    
    prelim = (targetFields * grandImaginaryMatrix * targetFields').*...,
        [integralMatrix_self,integralMatrix_mutual;integralMatrix_mutual,integralMatrix_self];
    
    poyntingMat = - R * prelim * R';
    
    %% Obtain flux from Poynting matrix
    poyntingMat = diag(real(poyntingMat(r1, r2)));
    for p = 1:N
        flux(p,i) = poyntingMat(p) + poyntingMat(N+p);
    end
    
    %% Calculate work input
    for j = 1:length(modLayer)
        ct = modLayer(j);
        S_matrix_mod = propagationSMatrix( thickness_list, M_matrix, f_matrix, numOfG, ct);
        [q_R, q_L] = meshgrid(diag(qValues{ct}));
        P_1 = (eye(2*N) - S_matrix{ct}(r1,r2) * S_matrix_mod{numOfLayer}(r2,r1)) \ S_matrix{ct}(r1,r1);
        P_2 = S_matrix_mod{numOfLayer}(r2,r1) * P_1;
        pad = [P_1 ; P_2] / Q_1 * [eye(2*N), Q_2];
        prelim_1 = pad * prelim * pad';
        
        if thickness_list(ct) == inf
            integralMatrix_self = 1./(1i*(q_L-conj(q_R)));
            integralMatrix_mutual = zeros(size(q_L));
        else
            integralMatrix_self = (1-exp(-1i*(q_L-conj(q_R))*thickness_list(ct)))./(1i*(q_L-conj(q_R)));
            integralMatrix_mutual = (exp(1i*conj(q_R)*thickness_list(ct))-exp(-1i*q_L*thickness_list(ct)))./(1i*(q_L+conj(q_R)));
        end
        
        prelim_2 = prelim_1.* [integralMatrix_self,integralMatrix_mutual;integralMatrix_mutual,integralMatrix_self];
        workMat = M_matrix{ct} * prelim_2 * M_matrix{ct}';
        
        workMat = workMat(r1,r1);
        
        epsilon = dielectricMatrix{ct};
        epsilon = epsilon - diag(diag(epsilon));
        work(i) = work(i) -imag(trace(kron(eye(2),W*epsilon) * workMat));
    end
    
end

%Normalize flux and work appropriately for if integrating over K
if strcmp(angle,'all')
    flux = flux/pi/2;
    work = work/pi/2;
end

end

