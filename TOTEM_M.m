function [] = TOTEM_M(filename)
    % Define the path to the source folder (change this as needed)
    srcFolder = 'src';
    addpath(srcFolder);
    % Create an instance of the InputReader class
    reader = InputReader(filename);
    fprintf('Initialized InputReader with filename: %s\n', inputfilename);
    mesh = Mesh(reader);
    bcinit = BCInit(reader, mesh);
end