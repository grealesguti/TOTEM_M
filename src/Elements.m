classdef Elements < handle
    properties
        % Define any properties needed for the class
    end
    
    methods
        % Constructor
        function obj = Elements()
            % Implement the constructor
        end
        
        % Destructor (MATLAB handles memory management, so no explicit destructor needed)
        
        % Evaluate shape functions for a linear triangular element with 4 nodes
        function EvaluateLinearTriangularShapeFunctions(obj, xi, eta, shapeFunctions)
            % Implement this method
        end
        
        function CalculateLinearTriangularShapeFunctionDerivatives(obj, xi, eta, shapeFunctiDerivatives)
            % Implement this method
        end

        % Evaluate shape functions for a linear quadrilateral element with 4 nodes
        function shapeFunctions = EvaluateLinearQuadrilateralShapeFunctions(~,xi, eta)
            % Define the shape functions as a 4x1 matrix
            shapeFunctions = zeros(4, 1);
        
            % Compute the shape functions for a linear quadrilateral element
            shapeFunctions(1) = 0.25 * (1 - xi) * (1 - eta); % N1
            shapeFunctions(2) = 0.25 * (1 + xi) * (1 - eta); % N2
            shapeFunctions(3) = 0.25 * (1 + xi) * (1 + eta); % N3
            shapeFunctions(4) = 0.25 * (1 - xi) * (1 + eta); % N4
        end
    
            
        function shapeFunctionDerivatives = EvaluateLinearQuadrilateralShapeFunctionDerivatives(~,xi, eta)
            % Ensure the shapeFunctionDerivatives matrix is of the correct size (4 nodes x 2 derivatives)
            shapeFunctionDerivatives = zeros(4, 2);
        
            % Define the derivatives of shape functions for a linear quadrilateral element
            shapeFunctionDerivatives(1, 1) = -0.25 * (1 - eta);  % dN1/dxi
            shapeFunctionDerivatives(1, 2) = -0.25 * (1 - xi);   % dN1/deta
        
            shapeFunctionDerivatives(2, 1) = 0.25 * (1 - eta);   % dN2/dxi
            shapeFunctionDerivatives(2, 2) = -0.25 * (1 + xi);  % dN2/deta
        
            shapeFunctionDerivatives(3, 1) = 0.25 * (1 + eta);   % dN3/dxi
            shapeFunctionDerivatives(3, 2) = 0.25 * (1 + xi);   % dN3/deta
        
            shapeFunctionDerivatives(4, 1) = -0.25 * (1 + eta);  % dN4/dxi
            shapeFunctionDerivatives(4, 2) = 0.25 * (1 - xi);   % dN4/deta
        end


        % Evaluate shape functions for a quadratic quadrilateral element with 8 nodes
        function shapeFunctions=EvaluateQuadraticQuadrilateralShapeFunctions(~, xi, eta)
            % Ensure the shapeFunctions vector is of the correct size (8 nodes)
            shapeFunctions = zeros(8, 1);
     
            % Shape functions for the four corner nodes
            shapeFunctions(1) = 0.25 * (1 - xi) * (1 - eta) * (- xi - eta - 1);
            shapeFunctions(2) = 0.25 * (1 + xi) * (1 - eta) * (+ xi - eta - 1);
            shapeFunctions(3) = 0.25 * (1 + xi) * (1 + eta) * (+ xi + eta - 1);
            shapeFunctions(4) = 0.25 * (1 - xi) * (1 + eta) * (- xi + eta - 1);
        
            % Shape functions for the four mid-side nodes
            shapeFunctions(5) = 0.5 * (1 - xi) * (1 + xi) * (1 - eta);
            shapeFunctions(6) = 0.5 * (1 + xi) * (1 + eta) * (1 - eta);
            shapeFunctions(7) = 0.5 * (1 - xi) * (1 + xi) * (1 + eta);
            shapeFunctions(8) = 0.5 * (1 - xi) * (1 - eta) * (1 - eta);
        end
        
        function shapeFunctionDerivatives=CalculateQuadraticQuadrilateralShapeFunctionDerivatives(~, xi, eta)
            % Ensure the shapeFunctionDerivatives matrix is of the correct size (8 nodes x 2 derivatives)
            shapeFunctionDerivatives = zeros(8, 2);
              
            shapeFunctionDerivatives(1, 1) = -0.25 * (-1 + eta) * (2*xi + eta);  % dN1/dxi
            shapeFunctionDerivatives(1, 2) = -0.25 * (-1 +xi) * (xi + 2*eta);   % dN1/deta

            shapeFunctionDerivatives(2, 1) = 0.25 * (-1 + eta) * (eta - 2*xi);   % dN2/dxi
            shapeFunctionDerivatives(2, 2) = 0.25 * (1 + xi) * (2*eta - xi);   % dN2/deta

            shapeFunctionDerivatives(3, 1) = 0.25 * (1 + eta) * (2*xi + eta);   % dN3/dxi
            shapeFunctionDerivatives(3, 2) = 0.25 * (1 + xi) * (xi + 2*eta);    % dN3/deta

            shapeFunctionDerivatives(4, 1) = -0.25 * (1 + eta) * (eta - 2*xi);  % dN4/dxi
            shapeFunctionDerivatives(4, 2) = -0.25 * (-1 + xi) * (2*eta -xi);    % dN4/deta

            shapeFunctionDerivatives(5, 1) = xi * (-1 + eta);  % dN5/dxi
            shapeFunctionDerivatives(5, 2) = 0.5 * (1 + xi) * (-1 + xi);    % dN5/deta

            shapeFunctionDerivatives(6, 1) = -0.5 * (1 + eta) * (-1 + eta);   % dN6/dxi
            shapeFunctionDerivatives(6, 2) = -eta * (1 + xi);    % dN6/deta

            shapeFunctionDerivatives(7, 1) = -xi * (1 + eta);    % dN7/dxi
            shapeFunctionDerivatives(7, 2) = -0.5 * (1 + xi) * (-1 - xi);   % dN7/deta
            
            shapeFunctionDerivatives(8, 1) = 0.5 * (1 + eta) * (-1 + eta);   % dN8/dxi
            shapeFunctionDerivatives(8, 2) = eta * (-1 + xi);   % dN8/deta
        end

        % Evaluate shape functions for a linear triangular element with 4 nodes
        function EvaluateLinearTetrahedraShapeFunctions(obj, xi, eta, zeta, shapeFunctions)
            % Implement this method
        end
        
        function CalculateLinearTetrahedraShapeFunctionDerivatives(obj, xi, eta, zeta, shapeFunctionDerivatives)
            % Implement this method
        end

        function nodes =  GetHexahedralNodeLocations(~)
            % Define the coordinates of the 8 nodes
            nodes = [
                -1, -1, -1;
                 1, -1, -1;
                 1,  1, -1;
                -1,  1, -1;
                -1, -1,  1;
                 1, -1,  1;
                 1,  1,  1;
                -1,  1,  1;
            ];

        end

        function nodes = GetSerendipityQuadHexahedralNodeLocations(~)
            % Define the coordinates of the 20 nodes
            nodes = [
                -1, -1, -1;
                 1, -1, -1;
                 1,  1, -1;
                -1,  1, -1;
                -1, -1,  1;
                 1, -1,  1;
                 1,  1,  1;
                -1,  1,  1;
                 0, -1, -1;
                 1,  0, -1;
                 0,  1, -1;
                -1,  0, -1;
                -1, -1,  0;
                 1, -1,  0;
                 1,  1,  0;
                -1,  1,  0;
                 0, -1,  1;
                 1,  0,  1;
                 0,  1,  1;
                -1,  0,  1;
            ];
        
        
        end        

        function shapeFunctions = EvaluateHexahedralLinearShapeFunctions(~,xi, eta, zeta)
            % Ensure the shapeFunctions vector is of the correct size (8 nodes)
            shapeFunctions = zeros(1, 8);
        
            % Define the shape functions for a linear hexahedral element
            shapeFunctions(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
            shapeFunctions(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
            shapeFunctions(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
            shapeFunctions(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
            shapeFunctions(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
            shapeFunctions(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
            shapeFunctions(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
            shapeFunctions(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
        end

        function shapeFunctionDerivatives = CalculateHexahedralLinearShapeFunctionDerivatives(~,xi, eta, zeta)
            % Ensure the shapeFunctionDerivatives matrix is of the correct size (8 nodes x 3 derivatives)
            shapeFunctionDerivatives = zeros(8, 3);
        
            % Define the derivatives of shape functions for a linear hexahedral element
            shapeFunctionDerivatives(1, 1) = -0.125 * (1 - eta) * (1 - zeta);  % dN1/dxi
            shapeFunctionDerivatives(1, 2) = -0.125 * (1 - xi) * (1 - zeta);   % dN1/deta
            shapeFunctionDerivatives(1, 3) = -0.125 * (1 - xi) * (1 - eta);    % dN1/dzeta
        
            shapeFunctionDerivatives(2, 1) = 0.125 * (1 - eta) * (1 - zeta);   % dN2/dxi
            shapeFunctionDerivatives(2, 2) = -0.125 * (1 + xi) * (1 - zeta);  % dN2/deta
            shapeFunctionDerivatives(2, 3) = -0.125 * (1 + xi) * (1 - eta);   % dN2/dzeta
        
            shapeFunctionDerivatives(3, 1) = 0.125 * (1 + eta) * (1 - zeta);   % dN3/dxi
            shapeFunctionDerivatives(3, 2) = 0.125 * (1 + xi) * (1 - zeta);   % dN3/deta
            shapeFunctionDerivatives(3, 3) = -0.125 * (1 + xi) * (1 + eta);  % dN3/dzeta
        
            shapeFunctionDerivatives(4, 1) = -0.125 * (1 + eta) * (1 - zeta);  % dN4/dxi
            shapeFunctionDerivatives(4, 2) = 0.125 * (1 - xi) * (1 - zeta);   % dN4/deta
            shapeFunctionDerivatives(4, 3) = -0.125 * (1 - xi) * (1 + eta);  % dN4/dzeta
        
            shapeFunctionDerivatives(5, 1) = -0.125 * (1 - eta) * (1 + zeta);  % dN5/dxi
            shapeFunctionDerivatives(5, 2) = -0.125 * (1 - xi) * (1 + zeta);   % dN5/deta
            shapeFunctionDerivatives(5, 3) = 0.125 * (1 - xi) * (1 - eta);    % dN5/dzeta
        
            shapeFunctionDerivatives(6, 1) = 0.125 * (1 - eta) * (1 + zeta);   % dN6/dxi
            shapeFunctionDerivatives(6, 2) = -0.125 * (1 + xi) * (1 + zeta);  % dN6/deta
            shapeFunctionDerivatives(6, 3) = 0.125 * (1 + xi) * (1 - eta);   % dN6/dzeta
        
            shapeFunctionDerivatives(7, 1) = 0.125 * (1 + eta) * (1 + zeta);   % dN7/dxi
            shapeFunctionDerivatives(7, 2) = 0.125 * (1 + xi) * (1 + zeta);   % dN7/deta
            shapeFunctionDerivatives(7, 3) = 0.125 * (1 + xi) * (1 + eta);   % dN7/dzeta
        
            shapeFunctionDerivatives(8, 1) = -0.125 * (1 + eta) * (1 + zeta);  % dN8/dxi
            shapeFunctionDerivatives(8, 2) = 0.125 * (1 - xi) * (1 + zeta);   % dN8/deta
            shapeFunctionDerivatives(8, 3) = 0.125 * (1 - xi) * (1 + eta);   % dN8/dzeta
        end

        % Evaluate shape functions for a hexahedral serendipity element with 20 nodes
        function shapeFunctions=CalculateHexahedralSerendipityShapeFunctions(~, xi, eta, zeta)
                % Ensure the shapeFunctions vector is of the correct size (20 nodes)
                shapeFunctions = zeros(1, 20);

                % Define the shape functions for a hexahedral serendipity element
                xi1 = xi;
                xi2 = eta;
                xi3 = zeta;
            
                shapeFunctions(1) = (1 - xi2) * (1 - xi3) * (1 - xi1) * (-xi2 - xi3 - xi1 - 2) / 8;
                shapeFunctions(2) = (1 + xi2) * (1 - xi3) * (1 - xi1) * (xi2 - xi3 - xi1 - 2) / 8;
                shapeFunctions(3) = (1 + xi2) * (1 + xi3) * (1 - xi1) * (xi2 + xi3 - xi1 - 2) / 8;
                shapeFunctions(4) = (1 - xi2) * (1 + xi3) * (1 - xi1) * (-xi2 + xi3 - xi1 - 2) / 8;
            
                shapeFunctions(5) = (1 - xi2) * (1 - xi3) * (1 + xi1) * (-xi2 - xi3 + xi1 - 2) / 8;
                shapeFunctions(6) = (1 + xi2) * (1 - xi3) * (1 + xi1) * (xi2 - xi3 + xi1 - 2) / 8;
                shapeFunctions(7) = (1 + xi2) * (1 + xi3) * (1 + xi1) * (xi2 + xi3 + xi1 - 2) / 8;
                shapeFunctions(8) = (1 - xi2) * (1 + xi3) * (1 + xi1) * (-xi2 + xi3 + xi1 - 2) / 8;
            
                shapeFunctions(9) = (1 - xi2^2) * (1 - xi3) * (1 - xi1) / 4;
                shapeFunctions(10) = (1 + xi2) * (1 - xi3^2) * (1 - xi1) / 4;
                shapeFunctions(11) = (1 - xi2^2) * (1 + xi3) * (1 - xi1) / 4;
                shapeFunctions(12) = (1 - xi2) * (1 - xi3^2) * (1 - xi1) / 4;
            
                shapeFunctions(13) = (1 - xi2^2) * (1 - xi3) * (1 + xi1) / 4;
                shapeFunctions(14) = (1 + xi2) * (1 - xi3^2) * (1 + xi1) / 4;
                shapeFunctions(15) = (1 - xi2^2) * (1 + xi3) * (1 + xi1) / 4;
                shapeFunctions(16) = (1 - xi2) * (1 - xi3^2) *(1 + xi1) / 4;
            
                shapeFunctions(17) = (1 - xi2) * (1 - xi3) * (1 - xi1^2) / 4;
                shapeFunctions(18) = (1 + xi2) * (1 - xi3) * (1 - xi1^2) / 4;
                shapeFunctions(19) = (1 + xi2) * (1 + xi3) * (1 - xi1^2) / 4;
                shapeFunctions(20) = (1 - xi2) * (1 + xi3) * (1 - xi1^2) / 4;

    
        end
        
        function shapeFunctionDerivatives=CalculateHexahedralSerendipityShapeFunctionDerivatives(~, xi, eta, zeta)
            xi1 = xi;
            xi2 = eta;
            xi3 = zeta;

            shapeFunctionDerivatives = zeros(3, 20);
        
            shapeFunctionDerivatives(1, 1) = ((xi2 - 1) * (xi3 - 1) * (xi1 + xi2 + xi3 + 2)) / 8 + ((xi1 - 1) * (xi2 - 1) * (xi3 - 1)) / 8;
            shapeFunctionDerivatives(1, 2) = -((xi2 + 1) * (xi3 - 1) * (xi1 - xi2 + xi3 + 2)) / 8 - ((xi1 - 1) * (xi2 + 1) * (xi3 - 1)) / 8;
            shapeFunctionDerivatives(1, 3) = ((xi2 + 1) * (xi3 + 1) * (xi1 - xi2 - xi3 + 2)) / 8 + ((xi1 - 1) * (xi2 + 1) * (xi3 + 1)) / 8;
            shapeFunctionDerivatives(1, 4) = -((xi2 - 1) * (xi3 + 1) * (xi1 + xi2 - xi3 + 2)) / 8 - ((xi1 - 1) * (xi2 - 1) * (xi3 + 1)) / 8;
            shapeFunctionDerivatives(1, 5) = ((xi1 + 1) * (xi2 - 1) * (xi3 - 1)) / 8 - ((xi2 - 1) * (xi3 - 1) * (xi2 - xi1 + xi3 + 2)) / 8;
            shapeFunctionDerivatives(1, 6) = -((xi2 + 1) * (xi3 - 1) * (xi1 + xi2 - xi3 - 2)) / 8 - ((xi1 + 1) * (xi2 + 1) * (xi3 - 1)) / 8;
            shapeFunctionDerivatives(1, 7) = ((xi2 + 1) * (xi3 + 1) * (xi1 + xi2 + xi3 - 2)) / 8 + ((xi1 + 1) * (xi2 + 1) * (xi3 + 1)) / 8;
            shapeFunctionDerivatives(1, 8) = -((xi2 - 1) * (xi3 + 1) * (xi1 - xi2 + xi3 - 2)) / 8 - ((xi1 + 1) * (xi2 - 1) * (xi3 + 1)) / 8;
            shapeFunctionDerivatives(1, 9) = -((xi2 * xi2 - 1) * (xi3 - 1)) / 4;
            shapeFunctionDerivatives(1, 10) = ((xi3 * xi3 - 1) * (xi2 + 1)) / 4;
            shapeFunctionDerivatives(1, 11) = ((xi2 * xi2 - 1) * (xi3 + 1)) / 4;
            shapeFunctionDerivatives(1, 12) = -((xi3 * xi3 - 1) * (xi2 - 1)) / 4;
            shapeFunctionDerivatives(1, 13) = ((xi2 * xi2 - 1) * (xi3 - 1)) / 4;
            shapeFunctionDerivatives(1, 14) = -((xi3 * xi3 - 1) * (xi2 + 1)) / 4;
            shapeFunctionDerivatives(1, 15) = -((xi2 * xi2 - 1) * (xi3 + 1)) / 4;
            shapeFunctionDerivatives(1, 16) = ((xi3 * xi3 - 1) * (xi2 - 1)) / 4;
            shapeFunctionDerivatives(1, 17) = -(xi1 * (xi2 - 1) * (xi3 - 1)) / 2;
            shapeFunctionDerivatives(1, 18) = (xi1 * (xi2 + 1) * (xi3 - 1)) / 2;
            shapeFunctionDerivatives(1, 19) = -(xi1*(xi2 + 1)*(xi3 + 1))/2;
            shapeFunctionDerivatives(1, 20) = (xi1*(xi2 - 1)*(xi3 + 1))/2;
        
            shapeFunctionDerivatives(2, 1) =((xi1 - 1)*(xi3 - 1)*(xi1 + xi2 + xi3 + 2))/8 + ((xi1 - 1)*(xi2 - 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(2, 2) =((xi1 - 1)*(xi2 + 1)*(xi3 - 1))/8 - ((xi1 - 1)*(xi3 - 1)*(xi1 - xi2 + xi3 + 2))/8;
            shapeFunctionDerivatives(2, 3) =((xi1 - 1)*(xi3 + 1)*(xi1 - xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(2, 4) =- ((xi1 - 1)*(xi3 + 1)*(xi1 + xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 - 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(2, 5) =- ((xi1 + 1)*(xi3 - 1)*(xi2 - xi1 + xi3 + 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(2, 6) =- ((xi1 + 1)*(xi3 - 1)*(xi1 + xi2 - xi3 - 2))/8 - ((xi1 + 1)*(xi2 + 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(2, 7) = ((xi1 + 1)*(xi3 + 1)*(xi1 + xi2 + xi3 - 2))/8 + ((xi1 + 1)*(xi2 + 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(2, 8) = ((xi1 + 1)*(xi2 - 1)*(xi3 + 1))/8 - ((xi1 + 1)*(xi3 + 1)*(xi1 - xi2 + xi3 - 2))/8;
            shapeFunctionDerivatives(2, 9) = -(xi2*(xi1 - 1)*(xi3 - 1))/2;
            shapeFunctionDerivatives(2, 10) = ((xi3*xi3 - 1)*(xi1 - 1))/4;
            shapeFunctionDerivatives(2, 11) = (xi2*(xi1 - 1)*(xi3 + 1))/2;
            shapeFunctionDerivatives(2, 12) = -((xi3*xi3 - 1)*(xi1 - 1))/4;
            shapeFunctionDerivatives(2, 13) = (xi2*(xi1 + 1)*(xi3 - 1))/2;
            shapeFunctionDerivatives(2, 14) = -((xi3^2 - 1)*(xi1 + 1))/4; %-((xi3^2 - 1)*(xi1 + 1))/4

            shapeFunctionDerivatives(2, 15) = -(xi2*(xi1 + 1)*(xi3 + 1))/2;
            shapeFunctionDerivatives(2, 16) = ((xi3*xi3 - 1)*(xi1 + 1))/4;
            shapeFunctionDerivatives(2, 17) = -((xi1*xi1 - 1)*(xi3 - 1))/4;
            shapeFunctionDerivatives(2, 18) = ((xi1*xi1 - 1)*(xi3 - 1))/4;
            shapeFunctionDerivatives(2, 19) = -((xi1*xi1 - 1)*(xi3 + 1))/4;
            shapeFunctionDerivatives(2, 20) = ((xi1*xi1 - 1)*(xi3 + 1))/4;
        
            shapeFunctionDerivatives(3, 1) =((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 + xi3 + 2))/8 + ((xi1 - 1)*(xi2 - 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(3, 2) =- ((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 + xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(3, 3) =((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 - xi3 + 2))/8 - ((xi1 - 1)*(xi2 + 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(3, 4) =((xi1 - 1)*(xi2 - 1)*(xi3 + 1))/8 - ((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 - xi3 + 2))/8;
            shapeFunctionDerivatives(3, 5) =- ((xi1 + 1)*(xi2 - 1)*(xi2 - xi1 + xi3 + 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 - 1))/8;
            shapeFunctionDerivatives(3, 6) =((xi1 + 1)*(xi2 + 1)*(xi3 - 1))/8 - ((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 - xi3 - 2))/8;
            shapeFunctionDerivatives(3, 7) =((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 + xi3 - 2))/8 + ((xi1 + 1)*(xi2 + 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(3, 8) =- ((xi1 + 1)*(xi2 - 1)*(xi1 - xi2 + xi3 - 2))/8 - ((xi1 + 1)*(xi2 - 1)*(xi3 + 1))/8;
            shapeFunctionDerivatives(3, 9) =-((xi2*xi2 - 1)*(xi1 - 1))/4;
            shapeFunctionDerivatives(3, 10) =(xi3*(xi1 - 1)*(xi2 + 1))/2;
            shapeFunctionDerivatives(3, 11) =((xi2*xi2 - 1)*(xi1 - 1))/4;
            shapeFunctionDerivatives(3, 12) =-(xi3*(xi1 - 1)*(xi2 - 1))/2;
            shapeFunctionDerivatives(3, 13) =((xi2*xi2 - 1)*(xi1 + 1))/4; 
            shapeFunctionDerivatives(3, 14) =-(xi3*(xi1 + 1)*(xi2 + 1))/2;
            shapeFunctionDerivatives(3, 15) =-((xi2*xi2 - 1)*(xi1 + 1))/4;
            shapeFunctionDerivatives(3, 16) =(xi3*(xi1 + 1)*(xi2 - 1))/2;
            shapeFunctionDerivatives(3, 17) =-((xi1*xi1 - 1)*(xi2 - 1))/4;
            shapeFunctionDerivatives(3, 18) =((xi1*xi1 - 1)*(xi2 + 1))/4;
            shapeFunctionDerivatives(3, 19) =-((xi1*xi1 - 1)*(xi2 + 1))/4;
            shapeFunctionDerivatives(3, 20) =((xi1*xi1 - 1)*(xi2 - 1))/4;

            shapeFunctionDerivatives=shapeFunctionDerivatives';
        end
        
    end
end
