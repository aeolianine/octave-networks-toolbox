% Find the "optimal" number of communities given a network using an eigenvector method
% Source: MEJ Newman: Finding community structure using the eigenvectors of matrices,
%                                                               arXiv:physics/0605087
% Newman, "Modularity and community structure in networks",
%                                                     arxiv.org/pdf/physics/0602124v1
%
% Q=(s^T)Bs, Bij=Aij-kikj/2m
% Bij^g=Bij - delta_ij * (sum k over g)B_ik
% Bij^g=(Aij-kikj/2m)-delta_ij(sum k over g)(A(g)_ik-deg(g)_i deg(k)_j/2(m_g))
% Bij^g=(Aij-kikj/2m)-delta_ij(k_i^(g)-k_i*sum(deg^(g)/2m)
%
% STEPS:
% 1 define current modularity matrix
% 2 compute eigenvector corresp. to largest eigenvalue
% 3 separate into 2 modules based on signs in eigenvector
% terminate when max eigenvalue is 0 for all subgraphs
%
% Other functions used: numEdges.m, degrees.m, subgraph.m, isConnected.m
% GB: last modified, Oct 12, 2012

function modules = newmanEigenvectorMethod(adj)

    modules = {};
    n = length(adj);
    m = numEdges(adj);
    [deg, ~, ~] = degrees(adj);
    queue{1} = [1:n]; % append all nodes to the queue

    while length(queue) > 0% while there is always a divisible subgraph

        % compute modularity matrix- Bg
        G = queue{1}; % nodes in current (sub)graph to partition
        adjG = subgraph(adj, G); % first adjG same as original adj

        [degG, ~, ~] = degrees(adjG);
        Bg = zeros(length(G));

        for i = 1:length(G) % nodes are G(i)

            for j = i:length(G) % nodes are G(j)

                % Bij = adj(G(i),G(j))-deg(G(i))*deg(G(j))/(2*m)
                % delta_ij = (i==j)
                % k_i^g - k_i (d_g)/2m = degG(i)-deg(G(i))*sum(degG)/(2*m)
                Bg(i, j) = (adj(G(i), G(j)) - deg(G(i)) * deg(G(j)) / (2 * m)) - (i == j) * (degG(i) - deg(G(i)) * sum(degG) / (2 * m));
                Bg(j, i) = Bg(i, j);

            end

        end

        [V, E] = eig(Bg); % eigenvalues

        if abs(max(diag(E))) <= 10^(-5)% terminate - indivisible
            queue = queue(2:length(queue));
            modules = {modules{1:length(modules)} G};
            continue
        end

        [~, indmax] = max(diag(E));
        u1 = V(:, indmax);
        comm1 = find(u1 > 0); comm2 = find(u1 < 0); % indices in G

        if not(isConnected(adj(G(comm1), G(comm1)))) || not(isConnected(adj(G(comm2), G(comm2))))

            queue = queue(2:length(queue)); % remove(G)
            modules = {modules{1:length(modules)} G}; % keep original G as indivisible
            % found disconnected subgraphs - do not explore this option

            continue

        end

        queue = {queue{1:length(queue)}, G(comm1), G(comm2)}; % add comm1 and comm2 to the queue
        queue = queue(2:length(queue));

    end



%!test
%!shared T
%! T = load_test_graphs();
%! modules = newmanEigenvectorMethod(T{4}{2});
%! assert(length(modules),2)
%! assert(modules{1},[4,5,6])
%! assert(modules{2},[1,2,3])

%!test
%! for x=1:50
%!   adj = randomGraph(randi(10)+5,1);
%!   Adj = zeros(4*length(adj));
%!   Adj(1:length(adj),1:length(adj))=adj;
%!   Adj(length(adj)+1:2*length(adj),length(adj)+1:2*length(adj))=adj;
%!   Adj(2*length(adj)+1:3*length(adj),2*length(adj)+1:3*length(adj))=adj;
%!   Adj(3*length(adj)+1:4*length(adj),3*length(adj)+1:4*length(adj))=adj;
%! 
%!   Adj(5,length(adj)+5)=1; Adj(length(adj)+5,5)=1; 
%!   Adj(length(adj)+6,2*length(adj)+6)=1; Adj(2*length(adj)+6,length(adj)+6)=1; 
%!   Adj(2*length(adj)+7,3*length(adj)+7)=1; Adj(3*length(adj)+7,2*length(adj)+7)=1; 
%!   Adj(3*length(adj)+1,1)=1; Adj(1,3*length(adj)+1)=1; 
%! 
%!   modules = newmanEigenvectorMethod(Adj);
%!   assert(length(modules),4)
%! 
%!   prescribed = randi(6)+2;
%!   
%!   n = randi(50)+50;
%!   adj = [];
%!   while not(isConnected(adj)); adj = randomModularGraph(n,prescribed,0.9*log(n)/n,2+0.3*rand); end
%!   modules = newmanEigenvectorMethod(adj);
%!   
%!   sumnodes = 0;
%!   for m=1:length(modules); sumnodes = sumnodes + length(modules{m}); end
%!   assert(sumnodes,n)
%!   
%!   for m1=1:length(modules)
%!     for m2=m1+1:length(modules)
%!       assert(length(intersect(modules{m1},modules{m2})),0)     
%!     end
%!   end
%!    
%! end


%!demo
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! modules = newmanEigenvectorMethod(bowtie)