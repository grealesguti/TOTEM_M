function qth = CalculateHeatFlow_nodes(reader,mesh,solver)
    total_number_of_nodes = length(mesh.data.NODE);
    qth=zeros(3,total_number_of_nodes);
    mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
    [node_el,etype,element_natural_coordinates] = mesh.retrievemeshtype(reader);
    initialdofs=solver.soldofs;
    for i=1:length(mesh_elements)
            elementTag= mesh_elements(i);
            xx= mesh.elements_density(elementTag);
            element_nodes = mesh.data.ELEMENTS{elementTag};
            number_of_nodes = length(element_nodes);
            element_coordinates=zeros(3,number_of_nodes);
            Tee=zeros(1,number_of_nodes);
            for nd=1:number_of_nodes
                Tee(nd)=initialdofs(element_nodes(nd)*2-1);
                element_coordinates(:,nd)=mesh.data.NODE{element_nodes(nd)};
            end

            element_material_index=mesh.elements_material(elementTag);

        for j=1:length(element_nodes)
                natural_coordinates=element_natural_coordinates(j,:);
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
                JM = dShape' * element_coordinates';
                DN = inv(JM) * dShape'; 
                Th = N * Tee';
                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                [Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                qth(:,element_nodes(j))=(-Dk*DN*Tee');
        end

    end
end

