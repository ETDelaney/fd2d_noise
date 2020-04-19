function [mask] = create_mask()
%CREATE_MASK Summary of this function goes here
%   create a tapered mask

    %- get configuration --------------------------------------------------
    [Lx, Lz, nx, nz] = input_parameters();
    [X, Z] = define_computational_domain(Lx, Lz, nx, nz);
    
    % load array file
    array_file = 'array_nref-172.mat';
    temp_array = load([fd2d_path(), 'output', filesep, array_file]);
    
    array_index = zeros([172 2]);
    
    % array is the 172 SWIM array configuration - convert to node
    array_index(:,1) = floor(temp_array.array(:,1)./(X(1,2)))+1;
    array_index(:,2) = floor(temp_array.array(:,2)./(Z(2,1)))+1;
    
    % now remove nodes in the middle of array
    array_index(47:127,:) = [];
    
    % populate new variables to define mask
    x = array_index(:,1);
    z = array_index(:,2);
    
    % create mask - notice the swap of z and x to match the gradient layout
    mask = poly2mask(z,x,125,141);
    
    % logical must be converted to double before applying a smoother
    mask = double(mask);
    
    % now use a Gaussian filter to make the mask smoother
    mask = imgaussfilt(mask,5);

end

