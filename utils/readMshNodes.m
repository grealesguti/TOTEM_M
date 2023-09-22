function nodes = readMshNodes(mshFileName)
    % Open the .msh file for reading
    fileID = fopen(mshFileName, 'r');
    
    if fileID == -1
        error('Failed to open the file.');
    end
    
    % Initialize an empty array to store nodes
    nodes = [];
    
    % Initialize a flag to identify the Nodes section
    inNodesSection = false;
    
    % Read lines from the file
    while ~feof(fileID)
        line = fgetl(fileID);
        
        % Check if we are in the Nodes section
        if strcmp(line, '$Nodes')
            inNodesSection = true;
            continue;
        end
        
        % Check if we are at the end of the Nodes section
        if strcmp(line, '$EndNodes')
            inNodesSection = false;
            break;
        end
        
        % If we are in the Nodes section, parse the node information
        if inNodesSection
            parts = strsplit(line);
            if numel(parts) == 4
                nodeID = str2double(parts{1});
                xCoord = str2double(parts{2});
                yCoord = str2double(parts{3});
                zCoord = str2double(parts{4});
                
                % Store the node information
                nodes = [nodes; [nodeID, xCoord, yCoord, zCoord]];
            end
        end
    end
    
    % Close the file
    fclose(fileID);
end
