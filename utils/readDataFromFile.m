function data = readDataFromFile(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    
    if fid == -1
        error('Could not open the file.');
    end
    
    data = struct();
    currentSection = '';
    sectionNames = {};  % Cell array to store section names
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
            if ~isfield(data, 'NODE')
                data.NODE = cell(1, 1);
            end
            
            data.NODE{nodeID} = nodeCoords;
            
            case '*ELEMENT'
                    % Read the element data
                    elementInfo = sscanf(line, '%d, %d, %d');
                    elementData = elementInfo(1:3);
                    elementtypes{numberofelements}=elementtype;
                    numberofelements=numberofelements+1;

                    % Store element data in the struct
                    if ~isfield(data, 'ELEMENTS')
                        data.ELEMENTS = {};
                    end
                    
                    data.ELEMENTS{end+1} = elementData;
            case '*ELSET'
                % Parse and store ELSET data
                elsetInfo = str2double(strsplit(line, ', '));
                elsetInfoWithoutNaN = elsetInfo(~isnan(elsetInfo));
                elsetArray = [elsetArray,elsetInfoWithoutNaN];
                % Store node data in the struct
                if ~isfield(data, 'ELSET')
                    data.ELSET = cell(1, 1);
                end
                
                data.ELSET{elsetID} = elsetArray;
                
            case '*NSET'  
                % Parse and store NSET data
                nsetInfo = str2double(strsplit(line, ', '));
                nsetInfoWithoutNaN = nsetInfo(~isnan(nsetInfo));
                nsetArray = [nsetArray,nsetInfoWithoutNaN];
                disp(nsetID);
                % Store node data in the struct
                if ~isfield(data, 'NSET')
                    data.NSET = cell(1, 1);
                end
                data.NSET{nsetID} = nsetArray;        
        end
    end

    data.ElementSelectionNames=elementselectionnames;
    data.NodalSelectionNames=nodalselectionnames;
    data.ElementTypes=elementtypes;
    % Close the file
    fclose(fid);
end
