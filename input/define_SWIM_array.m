function [array, platform] = define_SWIM_array(every4)

% this function allows us to build up the SWIM array for any computational
% domain - it also returns the position of the Oseberg C platform in
% reference to the array

load('src_rec_combos.mat')
[Lx, Lz, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_parameters();

domaincenter=[Lx/2 Lz/2];
arraycenter=[703.5 853.2]; % compute this with simpler trignometry - average of end points
centeringshift=domaincenter-arraycenter;


% Do you want every 4 receivers? Seems redundant to want more...
if strcmp(every4,'yes')
    %srclist=src_stations_every4;
    %reclist=rec_stations_every4;
    array=rec_stations_every4;
    
else
    %srclist=src_stations;
    %reclist=rec_stations;
    array=rec_stations;
    
end

% Center array to middle of domain
shiftingfactor=ones(1,2);
shiftingfactor=shiftingfactor.*centeringshift;
shiftingfactor=repmat(shiftingfactor,length(array(:,1)),1);

array(:,2:3) = array(:,2:3)+shiftingfactor;

platform = zeros(1,2);
platform(1,1) = array(1,2)+46.4; % it is 46.4 m east of rec 1
platform(1,2) = array(1,3)+145.45; % 145.45 m north of rec 1


