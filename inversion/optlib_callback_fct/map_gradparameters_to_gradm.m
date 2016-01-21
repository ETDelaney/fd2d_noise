
function [grad_m] = map_gradparameters_to_gradm(m, grad_parameters, usr_par)
% MAP_GRADPARAMETERS_TO_GRADM function to transform the partial derivatives
% w.r.t. the physical parameters to the partial derivatives
% w.r.t. to the model variable m. This requires to apply the chain rule
% which is implicitly given by the function MAP_M_TO_PARAMETERS.
%
% Input:
% m
% grad_parameters
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% grad_m
%
% See also MAP_M_TO_PARAMETERS.


if( strcmp( usr_par.type, 'source') )
    
        grad_m = reshape( imfilter( grad_parameters, usr_par.kernel.imfilter, 'circular' ), [], 1 );
    
elseif( strcmp( usr_par.type, 'structure') )
    
    % grad_m = 2 * usr_par.structure_inversion.v0^2 * reshape( usr_par.structure_inversion.rho .* grad_parameters, [], 1 ) .* (1+m);
    grad_m = 4.8e10 * reshape( imfilter( grad_parameters, usr_par.kernel.imfilter, 'circular' ), [], 1 );    
    
end


end