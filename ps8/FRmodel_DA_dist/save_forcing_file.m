% save_forcing_file

function save_forcing_file(forcing)

load control_params.mat
load forcing_original.mat % reload "full" forcing set

I = find(forcing(:,1)>=Day_beg & forcing(:,1)<=Day_end);
forcing = forcing(I(1)-3:I(end)+3,:);

% Save subset of interest
save forcing.mat forcing

return
