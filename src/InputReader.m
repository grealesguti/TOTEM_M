classdef InputReader < handle
    properties
        filename
        meshFileName
        MeshEntityName
        DesiredOutput
        physics
        MaterialProperties
        boundaryConditions
        elementMaterialIndices
        Material_index
        bctype
        bcloc
        bcval
        MaterialVolumes
        TopOpt_Objective
        TopOpt_ConstraintName
        TopOpt_CosntraintValue
        MMA_MaxIter
        TopOpt_ObjectiveSelection
        TopOpt_DesignElements
        TObctype
        TObcloc
        TObcminval
        TObcval
        TObcmaxval
    end
    
    methods
        function obj = InputReader(filename)
            obj.filename = filename;
            obj.TopOpt_Objective ='';
            obj.readFile();
            obj.addPenaltyKeys();
            %fprintf('Read Input: ' +filename+  '\n');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function readFile(obj)
            %fprintf('#INPUTREADER::readFile\n');
            inputFile = fopen(obj.filename, 'r');
            %disp(obj.filename)
            
            if inputFile == -1
                error('Error opening file.');
            end
            
            obj.meshFileName = 'defaultmesh.gmsh';
            
            while ~feof(inputFile)
                line = fgetl(inputFile);
                %%disp(line)
                tokens = strsplit(line);
                keyword = tokens{1};
                switch keyword
                    case 'mesh_file'
                        if numel(tokens) < 3
                            warning('Invalid mesh file and entity name.');
                        else
                            obj.meshFileName = tokens{2};
                            obj.MeshEntityName = tokens{3};
                            fprintf('New Mesh File and Mesh entity: %s %s\n', obj.meshFileName, obj.MeshEntityName);
                        end
                    case 'TopOpt_bc'
                            boundaryName = tokens{2};
                            surfaceName = tokens{3};
                            value = str2double(tokens{5});
                            minvalue = str2double(tokens{4});
                            maxvalue = str2double(tokens{6});
                            obj.bctype{end+1} =  boundaryName;
                            obj.bcloc{end+1} = surfaceName;
                            obj.bcval =[obj.bcval, value];
                            obj.TObctype{end+1} =  boundaryName;
                            obj.TObcloc{end+1} = surfaceName;
                            obj.TObcminval =[obj.TObcminval, minvalue];
                            obj.TObcval =[obj.TObcval, value];
                            obj.TObcmaxval =[obj.TObcmaxval, maxvalue];
                            fprintf('New TopOpt_bc entity: %s %s\n', boundaryName);
                    case 'MMA_MaxIter'
                            obj.MMA_MaxIter = tokens{2};
                            fprintf('New MMA_MaxIter: %s %s\n', obj.MMA_MaxIter);
                    case 'output'
                        if numel(tokens) < 2
                            warning('Invalid output keyword.');
                        else
                            obj.DesiredOutput = tokens{2};
                            fprintf('Output: %s\n', obj.DesiredOutput);
                        end
                    case 'physics'
                        if numel(tokens) < 2
                            warning('Invalid physics keyword.');
                        else
                            obj.physics = tokens{2};
                            fprintf('Physics: %s\n', obj.physics);
                        end
                    case 'material'
                        if numel(tokens) < 2
                            warning('Invalid material section.');
                        else
                            volumeName = tokens{2};
                            materialProperties = containers.Map();
                            
                            while ~feof(inputFile)
                                line = fgetl(inputFile);
                                if isempty(line) || strncmp(line, '##', 2)
                                    break;
                                end
                                
                                propertyTokens = strsplit(line);
                                if numel(propertyTokens) >= 2
                                    propertyName = propertyTokens{1};
                                    propertyPol = str2double(propertyTokens{2});
                                    propertyValue=zeros(1+propertyPol,1);
                                    for i = 3:(3+propertyPol)
                                        propertyValue(i-2)=str2double(propertyTokens{i});
                                    end
                                    materialProperties(propertyName) = propertyValue;
                                else
                                    warning('Invalid property format.');
                                    %fprintf('Line Content: %s\n', line);
                                end
                            end
                            
                            % Store material properties in MaterialProperties cell array
                            obj.MaterialProperties{end+1} = materialProperties;
                            obj.MaterialVolumes{end+1} = volumeName;
                        end
                    case 'bc'
                        if numel(tokens) < 4
                            warning('Invalid boundary condition format.');
                            %fprintf('Line Content: %s\n', line);
                        else
                            boundaryName = tokens{2};
                            surfaceName = tokens{3};
                            value = str2double(tokens{4});
                            obj.bctype{end+1} =  boundaryName;
                            obj.bcloc{end+1} = surfaceName;
                            obj.bcval =[obj.bcval, value];

                        end
                    case 'TopOpt_constraint'
                        if numel(tokens) < 4
                            warning('Invalid boundary condition format.');
                            %fprintf('Line Content: %s\n', line);
                        else
                            constraintName = tokens{2};
                            value = str2double(tokens{3});
                            obj.TopOpt_ConstraintName{end+1} =  constraintName;
                            obj.TopOpt_CosntraintValue =[obj.TopOpt_CosntraintValue, value];
                        end
                end
            end
            
            fclose(inputFile);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = getmaterialproperty(obj,index,propertyName)
            material = obj.MaterialProperties{index};
            if isKey(material, propertyName)
                value = material(propertyName);
            else
                disp(['Key "', propertyName, '" not found in the map.']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define a function to add Penalty keys
        function addPenaltyKeys(obj)
            % Loop through each materialProperties map
            for i = 1:numel(obj.MaterialProperties)
                material = obj.MaterialProperties{i};
                newMaterial = containers.Map();  % Create a new map for modified properties

                % Loop through each key-value pair in the current material
                keys = material.keys;
                for j = 1:numel(keys)
                    key = keys{j};

                    % Check if the key does not start with 'Penalty_'
                    if ~startsWith(key, 'Penalty_')
                        % Create a new key with 'Penalty_' prefix
                        newKey = ['Penalty_', key];

                        % Assign a value of 1 to the new key
                        newMaterial(newKey) = 1;
                    end
                end

                % Update the materialProperties map with the modified properties
                obj.MaterialProperties{i} = [material; newMaterial];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end
