function footprint = load_footprint_data(filename, varargin)

temp_data = load(filename);

field_name=fieldnames(temp_data);
spatial_footprints_cell = getfield(temp_data,field_name{1});

if nargin == 2
    cells2use = varargin{1};
else
    cells2use = 1:size(spatial_footprints_cell,2);
end

n_cells = length(cells2use);

if iscell(spatial_footprints_cell)
    footprint = zeros(n_cells,size(spatial_footprints_cell{1},1),size(spatial_footprints_cell{1},2));
    for cell_n = 1:n_cells
        cell_ID = cells2use(cell_n);
        footprint(cell_n,:,:) = spatial_footprints_cell{cell_ID};
    end
else
    footprint = spatial_footprints_cell;
end

end