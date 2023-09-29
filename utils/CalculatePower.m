function Power = CalculatePower(reader,mesh,solver)
    Power=0;
    integrationfunction = @(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) IntegrationPower(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
    mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
    for i=1:length(mesh_elements)
        element_tag= mesh_elements(i);
        [Element_Power,element_dof_indexes]=GeneralGaussIntegration(3, 14, element_tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_tag},integrationfunction);
        Power=Power+Element_Power;
    end
end

