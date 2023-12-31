function Power = CalculatePower(reader,mesh,solver)
    mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
    etype=mesh.data.ElementTypes{mesh_elements(1)};
    dim = mesh.retrieveelementdimension(etype); 
    integrationfunction = @(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx) IntegrationPower(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx);
    mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
    power_all=zeros(length(mesh_elements),1);
    parfor i=1:length(mesh_elements)
        element_tag= mesh_elements(i);
        [Element_Power,element_dof_indexes]=GeneralGaussIntegration(dim, reader.GI_order, element_tag, mesh, solver.soldofs,reader,mesh.data.ElementTypes{element_tag},integrationfunction);
        power_all(i)=Element_Power;
        %Power=Power+Element_Power;
    end
    Power=sum(power_all);
end

