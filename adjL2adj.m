% Convert an adjacency list to an adjacency matrix.
%
% INPUTS: adjacency list: length n, where L{i_1}=[j_1,j_2,...]
% OUTPUTS: adjacency matrix nxn
%
% Note: Assume that if node i has no neighbours, then L{i}=[];
% Last updated: May 20 2023

function adj = adjL2adj(adjL)

    adj = zeros(length(adjL));

    for i = 1:length(adjL)

        for j = 1:length(adjL{i})
            adj(i, adjL{i}(j)) = 1;
        end

    end

%!test
%!shared T
%! T = load_test_graphs();
%!assert(adjL2adj( T{9}{2} ),T{4}{2} )      % "bowtie" graph
%!assert(adjL2adj( T{17}{2} ),T{16}{2} )    % directed 3-cycle
%!assert(adjL2adj( T{12}{2} ), edgeL2adj(T{11}{2}) )


%!demo
%! aL = {[2, 3], [1, 3], [1, 2]};   % 3-cycle
%! adjL2adj(aL)
%! adjL2adj( { [2] ,[] } )   % single directed edge