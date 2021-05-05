classdef Mesh
    properties
        nNodes;
        nCells;
        nVar;
        nVertexElement;
        xGrid;
        yGrid;
        tri;
        vertexElements;
        neighbors;
        neighborsEdge;
        normals;
        edgeLengths;
        areas;
        xCenters;
        yCenters;
        incircles;
        data;
    end
    
    methods
        
        function this = Mesh(xL, xU, yL, yU, nX, nY, nVar)
            
            if nargin > 0
                
                % Create point cloud
                dx = (xU - xL) / nX;
                dy = (yU - yL) / nY;
                this.nNodes = 0;
                for i = 1:nX + 1
                    for j = 1:nY + 1
                        this.nNodes = this.nNodes + 1;
                        this.xGrid(this.nNodes) = xL + (i-1)*dx;
                        this.yGrid(this.nNodes) = yL + (j-1)*dy;
                    end
                end
                
                % Compute Delaunay triangulation
                this.tri = delaunay(this.xGrid, this.yGrid);
                [this.nCells, ~] = size(this.tri);
                
                % Compute the elements attached to each node
                this.nVertexElement = zeros(this.nNodes, 1);
                this.vertexElements = zeros(this.nNodes, 1);
                for i = 1:this.nCells
                    for k = 1:3
                        iNode = this.tri(i,k);
                        this.nVertexElement(iNode) = this.nVertexElement(iNode) + 1;
                        this.vertexElements(iNode, this.nVertexElement(iNode)) = i;
                    end
                end
                
                % Compute neighbors of each cell
                this.neighbors = zeros(this.nCells, 3);
                this.neighborsEdge = zeros(this.nCells, 3);
                edgeDef = [1 2; 2 3; 3 1];
                for i = 1:this.nCells
                    for iEdge = 1:3
                        iNode1 = this.tri(i, edgeDef(iEdge, 1));
                        iNode2 = this.tri(i, edgeDef(iEdge, 2));
                        
                        for k = 1:this.nVertexElement(iNode1)
                            j = this.vertexElements(iNode1, k);
                            if (i == j)
                                continue;
                            end
                            for jEdge = 1:3
                                jNode1 = this.tri(j, edgeDef(jEdge, 1));
                                jNode2 = this.tri(j, edgeDef(jEdge, 2));
                                if (iNode1 == jNode2 && iNode2 == jNode1)
                                    % Found common edge
                                    this.neighbors(i, iEdge) = j;
                                    this.neighborsEdge(i, iEdge) = jEdge;
                                    break;
                                end
                            end
                            if (this.neighbors(i, iEdge) > 0)
                                % Already found neighbor
                                break;
                            end
                        end
                    end
                end
                
                % Compute normal vectors, areas, and edge lengths
                z = [0; 0; 1];
                for i = 1:this.nCells
                    check = zeros(2,1); % aux variable to check mesh consistency
                    for iEdge = 1:3
                        iNode1 = this.tri(i,edgeDef(iEdge, 1));   % Global node number
                        iNode2 = this.tri(i,edgeDef(iEdge, 2));   % Global node number
                        v(:,iEdge) = [this.xGrid(iNode2) - this.xGrid(iNode1); this.yGrid(iNode2) - this.yGrid(iNode1); 0];
                        n3d = cross(v(:,iEdge),z);
                        this.normals(:,iEdge,i) = n3d(1:2) / sqrt(sum(n3d.^2));
                        this.edgeLengths(iEdge,i) = sqrt(sum(v(:,iEdge).^2));
                        check = check + this.normals(:,iEdge,i)*this.edgeLengths(iEdge,i);
                    end
                    if (sqrt(sum(check.^2)) > 1.e-12)
                        disp('Error: Mesh is not consistent!');
                        return;
                    end

                    % Compute area
                    temp = cross(v(:,1), -v(:,3));
                    this.areas(i) = temp(3) / 2;

                    % Barycenter of triangle
                    this.xCenters(i) = 1/3*(sum(this.xGrid(this.tri(i,:))));
                    this.yCenters(i) = 1/3*(sum(this.yGrid(this.tri(i,:))));

                    % Incircle diameter
                    this.incircles(i) = 4*this.areas(i)/sum(this.edgeLengths(:,i));
                end
                
                % Compute boundary cells
                for i = 1:this.nCells
                    for iEdge=1:3
                        iNode1 = this.tri(i, edgeDef(iEdge, 1));
                        iNode2 = this.tri(i, edgeDef(iEdge, 2));
                        
                        xNode1 = this.xGrid(iNode1);
                        xNode2 = this.xGrid(iNode2);
                        yNode1 = this.yGrid(iNode1);
                        yNode2 = this.yGrid(iNode2);
                        
                        if (xNode1 == xL && xNode2 == xL)
                            % West Boundary
                            this.neighbors(i, iEdge) = -1; 
                        elseif (xNode1 == xU && xNode2 == xU)
                            % East Boundary
                            this.neighbors(i, iEdge) = -2;
                        elseif (yNode1 == yL && yNode2 == yL)
                            % South Boundary
                            this.neighbors(i, iEdge) = -3;
                        elseif (yNode1 == yU && yNode2 == yU)
                            % North Boundary
                            this.neighbors(i, iEdge) = -4;
                        end
                    end
                end
                
                % Allocate space for data
                this.nVar = nVar;
                this.data = zeros(this.nVar, this.nCells);
                
            end
        end % End of Constructor
        
        function PlotMesh(this)
            figure;
            triplot(this.tri, this.xGrid, this.yGrid);
            title(sprintf('Mesh: nCells = %i, nNodes = %i', this.nCells, this.nNodes));
            xlabel('X');
            ylabel('Y');
        end
        
        function Plot(this)
            % Interpolate data from center to nodes
            nodeArea = zeros(this.nNodes, 1);
            nodeData = zeros(this.nVar, this.nNodes);
            for i = 1:this.nCells
                for k = 1:3
                    iNode = this.tri(i,k);
                    nodeArea(iNode) = nodeArea(iNode) + this.areas(i);
                    nodeData(:,iNode) = nodeData(:,iNode) + this.data(:,i)*this.areas(i);
                end
            end
            for iNode = 1:this.nNodes
                nodeData(:,iNode) = nodeData(:,iNode) / nodeArea(iNode);
            end
            
            for n = 1:this.nVar
                figure;
                s = trisurf(this.tri, this.xGrid, this.yGrid, nodeData(n,:));
                set(s, 'EdgeColor', [0 0 0]);
                set(s, 'FaceColor', 'interp');
                title(sprintf('Data Variable #%i', n));
                xlabel('X');
                ylabel('Y');
            end
        end
        
        function nodeData = CenterToNode(this)
            nodeArea = zeros(this.nNodes, 1);
            nodeData = zeros(this.nVar, this.nNodes);
            for i = 1:this.nCells
                for k = 1:3
                    iNode = this.tri(i,k);
                    nodeArea(iNode) = nodeArea(iNode) + this.areas(i);
                    nodeData(:,iNode) = nodeData(:,iNode) + this.data(:,i)*this.areas(i);
                end
            end
            for iNode = 1:this.nNodes
                nodeData(:,iNode) = nodeData(:,iNode) / nodeArea(iNode);
            end
        end
        
    end

end