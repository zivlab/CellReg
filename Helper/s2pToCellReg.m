function [footprint] = s2pToCellReg(s2p_path,varargin)
%S2PTOCELLREG 
% convert suite2p outputs to CellReg compatible footprint matrix 
% NoteI will later try to make the indeing faster (it is currently slow for
% what it is doing)

% Parameters:
%           - path of suite2p .mat file
%           - if pixels overlapping with surounding cells should be removed
%           from the footprint of the cell
% Output:
%           - footprint matrix (Ncells, ypix, xpix)

%% check inputs and load data

if nargin < 1 
    maskOverlap = varargin{1};
else 
    maskOverlap = 1;
end 

load(s2p_path,'stat','ops');

%% Process data 
nCells =length(stat);
footprint = zeros(nCells,ops.Lx, ops.Ly);

for itCell = 1:nCells
    % find subindices of the current neuron
    footprint_cell = zeros(size(footprint,2) , size(footprint,3) );
    idx = ...
        sub2ind([size(footprint,2) , size(footprint,3)] , ...
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
    footprint(itCell, :, :) = footprint_cell;
end 


end


