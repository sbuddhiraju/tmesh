function [ omega,epsilon ] = addMaterial( materialName )
%This function reads in the material name and output the material
%dielectric funcion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:    Kaifeng (Francis) Chen
% Copyright: Shanhui Fan's Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materialFileName=sprintf('%s%s',materialName,'.mat');
materialData=load(materialFileName);
omega=materialData.omega;
epsilon=materialData.epsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
