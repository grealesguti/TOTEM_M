classdef Mesh < handle
    properties 
        % Define class properties here
        % For example:
        inputReader
        data
        elements
        elements_material
        elements_density
        dofs_fixed
        element_characteristic_size
        Elements_volume
        Element_size
        dim
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
                    
                    %disp(line)        
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
                           %disp("Element Type:")
                           elementtype=tokens{2};
                           elementtype = strrep(elementtype, ' ', ''); % Remove spaces using strrep
                           elementtype=elementtype(6:end); % Remove "type="
                           %disp(elementtype)
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
                            %disp("node")
                        % Split the line and parse node data
                        nodeInfo = str2double(strsplit(line, ', '));
                        nodeID = nodeInfo(1);
                        nodeCoords = nodeInfo(2:end);
                        
                        % Store node data in the struct
                        if ~isfield(obj.data, 'NODE')
                            obj.data.NODE = cell(1, 1);
                        end
                            if strcmp(inputReader.Units,'mm')
                                nodeCoords=nodeCoords/1000;
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
                                elseif(elementtype=="T3D3")
                                    elementInfo = sscanf(line, '%d, %d, %d, %d');
                                    elementData = elementInfo(2:4);
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
                                elseif(elementtype=="CPS8")
                                    elementInfo = sscanf(line, '%d, %d, %d, %d, %d, %d, %d, %d, %d');
                                    elementData = elementInfo(2:9);
                                    elementIdx=elementInfo(1);
                                    elementtypes{elementIdx}=elementtype;
                                    numberofelements=numberofelements+1;        
                                elseif(elementtype=="C3D20")
                                    elementInfo = sscanf(line, '%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d,');
                                    line = fgetl(fid);
                                    elementInfo1 = sscanf(line, '%d, %d, %d, %d, %d');
                                    elementData(1:15) = elementInfo(2:16);
                                    elementData(16:20)= elementInfo1(1:5);
                                    elementIdx=elementInfo(1);
                                    elementtypes{elementIdx}=elementtype;
                                    numberofelements=numberofelements+1;     
                                else
                                    % Print a warning when an unknown element type is encountered
                                    warning('Unknown element type: %s', elementtype);
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
                            %disp(nsetID);
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
                obj.Elements_volume = obj.CalculateAllElementVolume(inputReader);

                if obj.dim==2
                    obj.Element_size=sqrt(min(obj.Elements_volume));
                elseif obj.dim==3
                    obj.Element_size=(min(obj.Elements_volume))^1/3;
                end
                % Close the file
                fclose(fid);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add getter and setter methods here
        function [shapeFunctions, shapeFunctionDerivatives] = selectShapeFunctionsAndDerivatives(~,etype, xi, eta, zeta)
            % Initialize output variables
            shapeFunctions = [];
            shapeFunctionDerivatives = [];
            elements=Elements();

            if etype == "CPS4" % 4-node quadrangle
                shapeFunctions = elements.EvaluateLinearQuadrilateralShapeFunctions(xi, eta);
                shapeFunctionDerivatives = elements.EvaluateLinearQuadrilateralShapeFunctionDerivatives(xi, eta);
            elseif etype == "CPS8" % 8-node second order quadrangle
                shapeFunctions = elements.EvaluateQuadraticQuadrilateralShapeFunctions(xi, eta);
                shapeFunctionDerivatives = elements.CalculateQuadraticQuadrilateralShapeFunctionDerivatives(xi, eta);
            elseif etype == "T3D2" % 8-node second order quadrangle
                shapeFunctions = elements.EvaluateLinearLineShapeFunctions(xi);
                shapeFunctionDerivatives = elements.EvaluateLinearLineShapeFunctionDerivatives(xi);
            elseif etype == "T3D3" % 8-node second order quadrangle
                shapeFunctions = elements.EvaluateQuadraticLineShapeFunctions(xi);
                shapeFunctionDerivatives = elements.EvaluateQuadraticLineShapeFunctionDerivatives(xi);
            elseif etype == "C3D8" % Hexahedral 8 node element
                shapeFunctions = elements.EvaluateHexahedralLinearShapeFunctions(xi, eta, zeta);
                shapeFunctionDerivatives = elements.CalculateHexahedralLinearShapeFunctionDerivatives(xi, eta, zeta);
            elseif etype == "C3D20" % Hexahedral 20 node element
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
            obj.elements_density = ones(1, length(obj.data.ELEMENTS));
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
        function elements = retrieveElementalSelection(obj,targetString)              
                % Use cellfun to find indices where targetString is located
                indices = find(cellfun(@(x) strcmp(x, targetString), obj.data.ElementSelectionNames));
                % Check if the length of indices is higher than 1 and issue a warning
                if length(indices) > 1
                    warning('Multiple occurrences of the target string found.');
                end
                elements = obj.data.ELSET{indices};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nodes = retrieveNodalSelection(obj,targetString)              

                % Use cellfun to find indices where targetString is located
                indices = find(cellfun(@(x) strcmp(x, targetString), obj.data.NodalSelectionNames));
                % Check if the length of indices is higher than 1 and issue a warning
                if length(indices) > 1
                    warning('Multiple occurrences of the target string found.');
                end
                nodes = obj.data.NSET{indices};
        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [node_el,etype_element,naturalcoordinates] = retrievemeshtype(obj,reader)
            elements_in_mesh = retrieveElementalSelection(obj, reader.MeshEntityName);
            prev_flag = 0;  % Initialize a previous flag
            cc=0;
            naturalcoordinates=[];
            %elements=Elements();
            for element = elements_in_mesh
                etype_element = obj.data.ElementTypes{element};
                if etype_element == "CPS4" % 4-node quadrangle
                    node_el = 4;
                    flag = 1;
                    naturalcoordinates=obj.elements.GetQuadNodeLocations();
                elseif etype_element == "CPS8" % 8-node second-order quadrangle
                    node_el = 8;
                    flag = 2;
                    naturalcoordinates=obj.elements.GetHexahedralNodeLocations();
                elseif etype_element == "C3D8" % Hexahedral 8 node element
                    node_el = 8;
                    flag = 3;
                    naturalcoordinates=obj.elements.GetHexahedralNodeLocations();
                elseif etype_element == "C3D20" % Hexahedral 20 node element
                    node_el = 20;
                    flag = 4;
                    naturalcoordinates=obj.elements.GetSerendipityQuadHexahedralNodeLocations();
                elseif etype_element == "T3D2" % 8-node second-order quadrangle
                    node_el = 2;
                    flag = 5;
                    naturalcoordinates=obj.elements.GetHexahedralNodeLocations();
                elseif etype_element == "T3D3" % 8-node second-order quadrangle
                    node_el = 3;
                    flag = 6;
                    naturalcoordinates=obj.elements.GetHexahedralNodeLocations();                    
                else
                    flag = 0;
                    % Handle unsupported element types or return an error code
                    % You can choose an appropriate error handling strategy here
                    % For example, you can throw an exception or set an error flag
                    % and handle it in the calling code.
                end
                
                % Check if the current flag is different from the previous flag
                if (cc>1 && flag ~= prev_flag)
                    warning('Different cases of etype_element found in the loop.');
                end
                
                % Update the previous flag
                prev_flag = flag;

                cc=cc+1;
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dim] = retrieveelementdimension(~,etype_element)


                if etype_element == "CPS4"      % DIM-2 4-node quadrangle
                    dim=2;      
                elseif etype_element == "T3D2"  % DIM-1 2-node line
                    dim=1;
                elseif etype_element == "T3D3"  % DIM-1 3-node line
                    dim=1;
                elseif etype_element == "CPS8"  % DIM-2 8-node second-order quadrangle
                    dim=2;
                elseif etype_element == "C3D8"  % DIM-3 8-node Hexahedral 8 node element
                    dim=3;
                elseif etype_element == "C3D20" % DIM-3 20-node Hexahedral 20 node element
                    dim=3;
                else
                    dim=-1;
                end
                          
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dim] = retrieveelementnumberofnodes(~,etype_element)


                if etype_element == "CPS4"      % DIM-2 4-node quadrangle
                    dim=4;      
                elseif etype_element == "T3D2"  % DIM-1 2-node line
                    dim=2;
                elseif etype_element == "T3D3"  % DIM-1 3-node line
                    dim=3;
                elseif etype_element == "CPS8"  % DIM-2 8-node second-order quadrangle
                    dim=8;
                elseif etype_element == "C3D8"  % DIM-3 8-node Hexahedral 8 node element
                    dim=8;
                elseif etype_element == "C3D20" % DIM-3 20-node Hexahedral 20 node element
                    dim=20;
                else
                    dim=-1;
                end
                          
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Elements_volume]=CalculateAllElementVolume(obj,reader)
            mesh_elements = obj.retrieveElementalSelection(reader.MeshEntityName);
                etype=obj.data.ElementTypes{mesh_elements(1)};
                obj.dim = obj.retrieveelementdimension(etype);             
            Elements_volume=zeros(length(mesh_elements),1);
            for i=1:length(mesh_elements)
                element_Tag=mesh_elements(i);
                element_nodes=obj.data.ELEMENTS{element_Tag};
                if obj.dim == 2
                    coordinates=zeros(3,4);
                    for j=1:4
                        coordinates(:,j)=obj.data.NODE{element_nodes(j)};
                    end                    
                    Elements_volume(i)=CalculateQuadArea(coordinates');
                else
                    coordinates=zeros(3,8);
                    for j=1:8
                        coordinates(:,j)=obj.data.NODE{element_nodes(j)};
                    end
                Elements_volume(i)=CalculateHexVolume(coordinates');
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
