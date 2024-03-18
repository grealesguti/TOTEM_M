function [N] = getTOEL(filepath)
        %filepath="TECTO/input_TECTO_Thermoel_Quad_cteF.txt";   
        reader = InputReader(filepath);

        mesh = Mesh(reader);
        TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
        N=length(TOEL);
end