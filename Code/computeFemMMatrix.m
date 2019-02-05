function [LHS] = computeFemMMatrix(xnodes, nconn)
ne = size(nconn,1);   % number of elements
nen = size(nconn,2);  % number of nodes per element
nn  = size(xnodes,1); % total number of nodes
nq  = 4;              % number of element integration points
ndim = 2;             % number of spatial dimensions (2 for 2D)

LHS = sparse(nn,nn);

% Loop over elements
for ielt = 1:ne
    % Element matrices
    LHSe = zeros(nen,nen);
    % Coordinates for element nodes
    coords = xnodes(nconn(ielt,:)',:);
    [~, wq, N, dNdx] = computeQuad2dFemShapeFunctions(coords); 
    % Form element matrix using gauss quadrature
    for iq = 1:nq
        % Loop over node pairs
        for i = 1:nen
            for j = 1:nen
                LHSe(i,j) =  LHSe(i,j)+N(iq,i)*N(iq,j)*wq(iq);  
            end
        end    
    end
    
    % Assemble to global matrix
    for i = 1:nen
        I = nconn(ielt,i);
        for j = 1:nen
            J = nconn(ielt,j);
            LHS(I,J) = LHS(I,J) + LHSe(i,j);
        end
    end
    
end
