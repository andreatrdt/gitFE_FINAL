function newLewis = gridExtract (oldGrid, oldLewis, newGrid)
% Extraction of the grid corresponding to the moneyness grid
%
%INPUT
% oldGrid:          grid vector
% oldLewis:         values vector
% newGrid:          grid in which the interpolation is needed

    newLewis = interp1(oldGrid, oldLewis, newGrid);

end % function gridExtract