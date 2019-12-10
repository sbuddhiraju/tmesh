function [flux, work] = integrateK(  omega, thickness_list,dielectricMatrix, dielectricMatrixInverse, dielectric_im_Matrix, numOfG, sourceLayer,targetLayer,modLayer, period,angle)
%This function will do the integral over ky
% kx_list                                    the parellel kx vector, normalized
% omega:                                     angular frequency
% dielectricMatrix, dielectricMatrixInverse: the matrices used to calculate
%                                            the eigenmodes
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

if strcmp(angle, 'all')
    numOfPoints = 400;
    upperBound = 2; %>1 to accommodate higher sidebands
    k_list = (0:upperBound/numOfPoints:upperBound);

else
    k_list = sin(pi/180 * 0);
end

[flux,work] = fluxAtKxKy( omega, thickness_list, 0, k_list,dielectricMatrix,...
    dielectricMatrixInverse, dielectric_im_Matrix, numOfG, sourceLayer,targetLayer,modLayer, period,angle);


if strcmp(angle,'all')
    flux = sum(flux.*k_list*omega,2) * (omega*upperBound/numOfPoints); %the omega here is just a scale factor for k_par, has no frequency significance.
    work = sum(work.*k_list*omega,2) * (omega*upperBound/numOfPoints); %the omega here is just a scale factor for k_par, has no frequency significance.
    
end

flux =  2 * omega * flux / pi; %this is from <JJ*>
work = 2 * omega * work / pi; %this is from <JJ*>
%the omega here is the frequency of the thermal source.

end

