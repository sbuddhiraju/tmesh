%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes thermal emission into vacuum from planar, layered
% media whose dielectric functions are time-dependent.

% tMESH: Siddharth Buddhiraju
% Based on MESH by Kaifeng Chen
% Copyright: Shanhui Fan's Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
PhysicsConst;

%% Preparing data files
omega = linspace(156,186,500)'*1e+12;

epsilon = (14-1e-10i)*ones(length(omega),1);
save Dielectric1.mat epsilon omega;
epsilon = (4-1e-10i)*ones(length(omega),1);
save Dielectric2.mat epsilon omega;

ep_inf = 6.7;
w_lo = 1.83e+14;
w_to = 1.49e+14;
gamma = 8.97e+11;
frac = 0.1;
ep2 = 14;
epsilon = ep_inf * (w_lo^2-omega.^2+1i*omega*gamma)./(w_to^2-omega.^2+1i*omega*gamma);
epsilon = frac*epsilon + (1-frac)*ep2;
save SiC.mat epsilon omega;

epsilon = (1-1e-10i)*ones(length(omega),1);
save Vacuum.mat epsilon omega;
epsilon = (1-1i)*1e+10*ones(length(omega),1);
save PEC.mat epsilon omega;

%% Define structure
materialNames={'Vacuum';'Dielectric1';'Dielectric2';'PEC';'SiC'};
%list of materials used

numOfG = 1; %number of sidebands

bragg1 = [2,0;3,0]; layers1 = 6; 
bragg2 = [3,1;2,0]; layers2 = 3;
bragg3 = [3,0;2,0]; layers3 = 3;    
%pattern and number of repetitions of Bragg layers

sysInfo = [4,0; repmat(bragg1,layers1,1);
    5,0;repmat(bragg2,layers2,1);
    3,1;2,1;
    repmat(bragg3,layers3,1);
    1,0];
%first element is material number from materilNames, 
%second is modulation (1 or 0)

strength = 0.1; %modulation strength
Omega = 1.0381e+13; %modulation frequency
period = 2 * pi / Omega; %temporal period

scale = 1e-6;
bragg_thickness = [1;1]; %thickness of Bragg layers
def1 = 2.1; %thickness of Defect 1
def2 = 1.6; %thickness of Defect 2

thickness_list=scale * [Inf;repmat(bragg_thickness,layers1,1);def1;...
    repmat(bragg_thickness,layers2,1);bragg_thickness(2);def2;...
    repmat(bragg_thickness,layers3,1);Inf];


%% Set related variables
numOfLayer=length(thickness_list);
numOfMaterials=length(materialNames);
sourceLayer=1 + 2*layers1 + 1; %index of source layer  
%make sure the source layer is a lossy layer
modLayer = find(sysInfo(:,2)==1); %this is where work is computed
targetLayer = length(thickness_list)-1; %flux is computed here
%targetLayer should be lossless

angle = 'normal'; % 'normal' for single channel or 'all' for integration
console = 'on'; %progress bar on console output
epsilon=cell(numOfMaterials,1);
N = 2*numOfG + 1;

for i=1:numOfMaterials
    [omega_list,epsilon{i}]=addMaterial(materialNames{i});
end

omega_list = omega_list/c;

%% Begin computation
fluxSpectrum = zeros(N,length(omega_list));
work = zeros(1,length(omega_list));

for i = 1:length(omega_list)
    if mod(i,10) == 0 && strcmp(console,'on')
        percentDone = 100 * i/length(omega_list);
        msg = sprintf('Computing emission spectrum: %3.1f', percentDone); %Don't forget this semicolon
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
    end
    
    epsilonMatrix = cell(numOfLayer,1);
    epsilonMatrixInverse = cell(numOfLayer,1);
    epsilon_im_Matrix = cell(numOfLayer,1);
    omega = omega_list(i);
    
    for j = 1:numOfLayer
        index = sysInfo(j,1);
        if sysInfo(j,2) == 0
            epsilonMatrix{j} = epsilon{index}(i) * eye(N);
            epsilonMatrixInverse{j} = 1/epsilon{index}(i) * eye(N);
            epsilon_im_Matrix{j} = imag(epsilon{index}(i)) * eye(N);
            
        elseif sysInfo(j,2) == 1
            epsilonMatrix{j} = epsilon{index}(i) * eye(N) + strength / 2 * ...
                (diag(ones(1,N-1),1) + diag(ones(1,N-1),-1));
            epsilonMatrixInverse{j} = inv(epsilonMatrix{j});
            epsilon_im_Matrix{j} = imag(epsilon{index}(i)) * eye(N);
        else
            disp('Invalid modulation status in sysInfo. Exit.'); return;
            
        end
    end
    
    %% Compute heat flux
    [fluxSpectrum(:,i), work(i)] = integrateK(  omega, thickness_list, epsilonMatrix, epsilonMatrixInverse, epsilon_im_Matrix, numOfG, sourceLayer,targetLayer,modLayer,period,angle);
    
end

%% Saving data
omega_list = omega_list*c;
filename = 'output.mat';
save(filename);

%% [Optional] Process and plot saved data
load(filename);

T = 290; %Set temperature to create reference blackbody spectrum
%This is useful to compute emissivity only in the presence of modulation
TH = @(omega,T) hbar*omega./(exp(hbar*omega/kB/T)-1); %thermal distribution

% Defining book-keeping variables
fS = fluxSpectrum;
domega = omega_list(2)-omega_list(1);
shift = round(Omega/domega);
flux = zeros(1,2*numOfG*shift + length(omega_list));
fS = fS.*TH(omega_list,T)';
for n = -numOfG:numOfG
    start = 1 + (n+numOfG)*shift;
    flux(start:start+length(omega_list)-1) = ...
        flux(start:start+length(omega_list)-1) + fS(n + numOfG + 1,:);
end
omega_shift = numOfG*shift*domega;
omega_list_new = linspace(omega_list(1)-omega_shift,omega_list(end)+...
    omega_shift,length(flux));
domega_new = omega_list_new(2)-omega_list_new(1);

%Define density of states
if strcmp(angle,'all')
    factor = omega_list_new.^2/4/pi^2/c^2; %all propagating channels
else
    factor = 1/pi; %single propagating channel
end

%Define emissivity
emissivity = flux./TH(omega_list_new,T)/factor;

%Plot
figure;
plot(omega_list_new/(2*pi*1e+12),emissivity,'linewidth',2); hold on;
set(gca,'linewidth',1,'fontsize',15);
xlabel('Frequency (THz)');
ylabel('\epsilon(\omega, k = 0)');