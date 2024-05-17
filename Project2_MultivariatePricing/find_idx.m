function idx = find_idx(data, date)
% Computation that finds the index in a struct given a certain date
%
% INPUT:
% data: struct containing the data
% date: date to be found
% 
% OUTPUT:
% idx: index of the date in the struct

    idx = find(data == date);

end % function find_idx