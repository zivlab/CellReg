classdef get_spatial_footprints
    properties
        footprints
        path
        size
        write2path
    end
    methods 
        function obj = get_spatial_footprints(val)
            obj.path = val;
            obj = get_size(obj);
            obj = get_file_type(obj);
        end
        function obj = load_footprints(obj, varargin)
            if isempty(obj.footprints)
                if strmatch(class(obj.path),'char')
                    if nargin == 2
                        cell2use = varargin{1};
                        obj.footprints = load_footprint_data(obj.path, cell2use);
                    else
                        obj.footprints = load_footprint_data(obj.path);
                    end
                else
                    obj.footprints = obj.path;
                    obj = remove_path(obj);
                end
            end
        end
        function obj = remove_path(obj)
            obj.path = [];
        end
        function obj = get_size(obj)
            if strmatch(class(obj.path),'char')
                matObj = matfile(obj.path);
                deets = whos(matObj);
                if length(deets.size) < 3
                    ncells = eval('size(matObj.(deets.name),2)');
                    dims_FOV = eval('matObj.(deets.name)(1,1)');
                    dims_FOV = size(dims_FOV{1});
                    obj.size = [ncells, dims_FOV];
                else
                    obj.size = deets.size;
                end
            else
                obj.size = size(obj.path);
            end
        end
        function obj = get_file_type(obj)
            if strmatch(class(obj.path),'char')
                temp_save_loc = strsplit(obj.path,filesep);
                save_loc = temp_save_loc{1};
                for part = 2:length(temp_save_loc) -1
                    save_loc = [save_loc, filesep, temp_save_loc{part}];
                end
                obj.write2path = save_loc;
            else
                obj.write2path = [];
            end
        end
       
    end
end