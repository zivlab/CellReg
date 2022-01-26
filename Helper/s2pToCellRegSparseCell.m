function [footprint] = s2pToCellRegSparseCell(s2p_path,varargin)
%S2PTOCELLREG 
% convert suite2p outputs to CellReg compatible footprint cell containing
% sparse footprints for each cell

% Parameters:
%           - path of suite2p .mat file or the file itself
%           - if pixels overlapping with surounding cells should be removed
%           from the footprint of the cell
% Output:
%           - footprint cell array (Ncells x 1) where each cell element
%           contains sparse matrix (ypix x xpix) containing footprint data
%           for that ROI 

%% check inputs and load data

if nargin < 1 
    maskOverlap = varargin{1};
else 
    maskOverlap = 1;
end 

if ~isstruct(s2p_path)
    load(s2p_path,'stat','ops');
    nCells =length(stat);
    cell_ns = 1:nCells;
else
    stat = s2p_path.stat;
    ops = s2p_path.ops;
    nCells = length(s2p_path.iscell);
    cell_ns = s2p_path.iscell;
end
%% Process data 
% footprint = zeros(nCells,ops.Lx, ops.Ly);

for cell_n = 1:nCells
    itCell = cell_ns(cell_n);
    % find subindices of the current neuron
    footprint_cell = zeros(ops.Lx , ops.Ly );
    idx = ...
        sub2ind([size(footprint_cell,1) , size(footprint_cell,2)] , ...
        stat{1,itCell}.ypix,stat{1,itCell}.xpix );
    
    % remove pixel that overlap with surounding cells
    if maskOverlap
        idx = idx(~stat{1,itCell}.overlap);
        lam = stat{1,itCell}.lam(~stat{1,itCell}.overlap);
    else 
        lam = stat{1,itCell}.lam;
    end
    
    % update the footprint matrix
    footprint_cell(idx) =  lam;
    footprint{cell_n} = sparse(footprint_cell);
end 


end
