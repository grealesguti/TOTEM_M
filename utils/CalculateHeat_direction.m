function Q_dir = CalculateHeat_direction(reader,mesh,solver,surface,direction)
    Q=zeros(3,1);
    integrationfunction = @(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) IntegrationHeat(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
    mesh_elements = mesh.retrieveElementalSelection(surface);
    for i=1:length(mesh_elements)
        element_tag= mesh_elements(i);
        [Element_Heat,element_dof_indexes]=GeneralGaussIntegration(2, 5, element_tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_tag},integrationfunction);
        Q=Q+Element_Heat;
    end
    Q_dir=Q*direction;
end

