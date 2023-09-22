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
    end
    
    methods
        function obj = InputReader(filename)
            obj.filename = filename;
            obj.readFile();
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
                                    propertyValue = str2double(propertyTokens{2});
                                    materialProperties(propertyName) = propertyValue;
                                else
                                    warning('Invalid property format.');
                                    %fprintf('Line Content: %s\n', line);
                                end
                            end
                            
                            % Store material properties in MaterialProperties cell array
                            obj.MaterialProperties{end+1} = materialProperties;
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
                end
            end
            
            fclose(inputFile);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end
