classdef Mesh < handle
    properties 
        % Define class properties here
        % For example:
        inputReader
        data
        elements
        elements_material
        dofs_fixed
    end
    
    methods
        % Constructor
        function obj = Mesh(inputReader)
            % Constructor code here
            % Initialize class properties based on input
            obj.elements=Elements();
            obj.ReadMesh(inputReader);    
            obj.DefineMaterials(inputReader)
        end
        
        % Add other methods here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ReadMesh(obj,inputReader)
            % Open the .msh file for reading
                % Open the file for reading
                fid = fopen(inputReader.meshFileName, 'r');
                
                if fid == -1
                    error('Could not open the file.');
                end
                
                obj.data = struct();
                currentSection = '';
                elementtypes = {};    % Loop through the lines in the file
                elementselectionnames = {};
                nodalselectionnames = {};
                numberofelements=1;
                nsetID=0;
                elsetID=0;
               while ~feof(fid)
                    line = fgetl(fid);
                    
                    disp(line)        
                    % Check if the line is a section heading and if there are commas
                    if startsWith(line, '*')
                        if contains(line, ',')
                            % If there are commas, split it at the comma
                            tokens = strsplit(line, ',');
                            %disp(tokens)
                            currentSection = tokens{1};
                        else
                            currentSection=line;
                        end
                        disp("Current Section:")
                        disp(currentSection)
                        % Extract and store the section name
                        if (currentSection=="*ELEMENT")
                           disp("Element Type:")
                           elementtype=tokens{2};
                           elementtype = strrep(elementtype, ' ', ''); % Remove spaces using strrep
                           elementtype=elementtype(6:end); % Remove "type="
                           disp(elementtype)
                        elseif (currentSection=="*ELSET")
                           disp("Selection Name:")
                           selectionname=tokens{2};
                           selectionname=selectionname(7:end); % Remove "ELSET="
                           disp(selectionname)  
                           elsetID=elsetID+1;
                           elsetArray = []; % Example existing array
                           elementselectionnames{elsetID}=selectionname;
                        elseif (currentSection=="*NSET")
                           disp("Selection Name:")
                           selectionname=tokens{2};
                           selectionname=selectionname(6:end); % Remove "NSET="
                           disp(selectionname)      
                           nsetID = nsetID+1;
                           nsetArray = []; % Example existing array
                           nodalselectionnames{nsetID}=selectionname;
                        end
                        continue
                    end
                    % Parse data based on the current section
                    switch currentSection
                    case '*NODE'
                            disp("node")
                        % Split the line and parse node data
                        nodeInfo = str2double(strsplit(line, ', '));
                        nodeID = nodeInfo(1);
                        nodeCoords = nodeInfo(2:end);
                        
                        % Store node data in the struct
                        if ~isfield(obj.data, 'NODE')
                            obj.data.NODE = cell(1, 1);
                        end
                        
                        obj.data.NODE{nodeID} = nodeCoords;
                        
                        case '*ELEMENT'
                                % Read the element data
                                if(elementtype=="T3D2") 
                                    elementInfo = sscanf(line, '%d, %d, %d');
                                    elementData = elementInfo(2:3);
                                    elementIdx=elementInfo(1);
                                    elementtypes{elementIdx}=elementtype;
                                    numberofelements=numberofelements+1;
                                elseif(elementtype=="CPS4")
                                    elementInfo = sscanf(line, '%d, %d, %d, %d, %d');
                                    elementData = elementInfo(2:5);
                                    elementIdx=elementInfo(1);
                                    elementtypes{elementIdx}=elementtype;
                                    numberofelements=numberofelements+1;
                                elseif(elementtype=="C3D8")
                                    elementInfo = sscanf(line, '%d, %d, %d, %d, %d, %d, %d, %d, %d');
                                    elementData = elementInfo(2:9);
                                    elementIdx=elementInfo(1);
                                    elementtypes{elementIdx}=elementtype;
                                    numberofelements=numberofelements+1;
                                else
                                end
            
                                % Store element data in the struct
                                if ~isfield(obj.data, 'ELEMENTS')
                                    obj.data.ELEMENTS = {};
                                end
                                
                                obj.data.ELEMENTS{elementIdx} = elementData;
                        case '*ELSET'
                            % Parse and store ELSET data
                            elsetInfo = str2double(strsplit(line, ', '));
                            elsetInfoWithoutNaN = elsetInfo(~isnan(elsetInfo));
                            elsetArray = [elsetArray,elsetInfoWithoutNaN];
                            % Store node data in the struct
                            if ~isfield(obj.data, 'ELSET')
                                obj.data.ELSET = cell(1, 1);
                            end
                            
                            obj.data.ELSET{elsetID} = elsetArray;
                            
                        case '*NSET'  
                            % Parse and store NSET data
                            nsetInfo = str2double(strsplit(line, ', '));
                            nsetInfoWithoutNaN = nsetInfo(~isnan(nsetInfo));
                            nsetArray = [nsetArray,nsetInfoWithoutNaN];
                            disp(nsetID);
                            % Store node data in the struct
                            if ~isfield(obj.data, 'NSET')
                                obj.data.NSET = cell(1, 1);
                            end
                            obj.data.NSET{nsetID} = nsetArray;        
                    end
                end
            
                obj.data.ElementSelectionNames=elementselectionnames;
                obj.data.NodalSelectionNames=nodalselectionnames;
                obj.data.ElementTypes=elementtypes;
                % Close the file
                fclose(fid);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add getter and setter methods here
        function [shapeFunctions, shapeFunctionDerivatives] = selectShapeFunctionsAndDerivatives(~,etype, xi, eta, zeta)
            % Initialize output variables
            elements=Elements();
            shapeFunctions = [];
            shapeFunctionDerivatives = [];
        
            if etype == "CPS4" % 4-node quadrangle
                shapeFunctions = elements.EvaluateLinearQuadrilateralShapeFunctions(xi, eta);
                shapeFunctionDerivatives = elements.EvaluateLinearQuadrilateralShapeFunctionDerivatives(xi, eta);
            elseif etype == 16 % 8-node second order quadrangle
                shapeFunctions = elements.EvaluateQuadraticQuadrilateralShapeFunctions(xi, eta);
                shapeFunctionDerivatives = elements.CalculateQuadraticQuadrilateralShapeFunctionDerivatives(xi, eta);
            elseif etype == "C3D8" % Hexahedral 8 node element
                shapeFunctions = elements.EvaluateHexahedralLinearShapeFunctions(xi, eta, zeta);
                shapeFunctionDerivatives = obj.elements.CalculateHexahedralLinearShapeFunctionDerivatives(xi, eta, zeta);
            elseif etype == 17 % Hexahedral 20 node element
                shapeFunctions = elements.CalculateHexahedralSerendipityShapeFunctions(xi, eta, zeta);
                shapeFunctionDerivatives = elements.CalculateHexahedralSerendipityShapeFunctionDerivatives(xi, eta, zeta);
            else
                % Handle unsupported element types or return an error code
                % You can choose an appropriate error handling strategy here
                % For example, you can throw an exception or set an error flag
                % and handle it in the calling code.
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Example getter:
        function DefineMaterials(obj,inputReader)
            % Initialize the output
            obj.elements_material = zeros(1, length(obj.data.ELEMENTS));
        
            % Loop through each material property
            for i = 1:length(inputReader.MaterialProperties)
                % Get the ElementSelectionName for this material
                elementSelectionName = inputReader.MaterialVolumes{i};
        
                % Find the corresponding ELSET component
                flag=0;
                if isfield(obj.data, 'ELSET') && isfield(obj.data, 'ElementSelectionNames')
                    elsetNames = obj.data.ElementSelectionNames;
                    for j = 1:length(elsetNames)
                        if strcmp(elementSelectionName, elsetNames{j})
                            obj.elements_material(obj.data.ELSET{j}) = i;
                            flag=1;
                            break;
                        end
                    end
                else
                    % Display a warning if ELSET or ElementSelectionNames are missing
                    warning('ELSET or ElementSelectionNames data not found.');
                end
        
                % Display a warning if no matching ELSET component is found
                if flag == 0
                    warning(['No matching ELSET component found for material ' ...
                        'name %d (%s).'], i, elementSelectionName);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end
