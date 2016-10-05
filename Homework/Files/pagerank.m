function x = pagerank(U,G,p)
% PAGERANK  Google's PageRank
% x = pagerank(U,G,p) uses the URLs and adjacency matrix produced by SURFER,
% together with a damping factory p, (default is .85), to compute the
% page rank, and print the dominant URLs in page rank order.

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree, r = in-degree

[n,n] = size(G);
c = sum(G,1);

% Any node with no out links is temporarily linked to all nodes

for j = find(c==0);
   G(:,j) = 1;
end
   
% Solve (I - p*G*D)x = delta*e.

if nargin < 3, p = .85; end
delta = (1-p)/n;
e = ones(n,1);
D = spdiags(1./sum(G)',0,n,n);
I = speye(n,n);
x = (I - p*G*D)\(delta*e);

% Print the top 10 URLs in page rank order.
[~,id] = sort(x,1,'descend');
fprintf('#     PageRank     Page\n'); 
for j=1:10
   fprintf('%02d    %1.2e     %s\n',j,x(id(j)),U{id(j)});
end

[~,id] = sort(x,1,'ascend');
fprintf('\n#     PageRank     Page\n'); 
for j=1:10
   fprintf('%02d    %1.2e     %s\n',n - 10 + j,x(id(j)),U{id(j)});
end

end