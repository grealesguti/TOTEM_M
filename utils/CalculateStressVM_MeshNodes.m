function [sigmaVM0_ordered] = CalculateStressVM_MeshNodes(reader,mesh,solver)

mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
total_number_of_elements = length(mesh_elements);
total_number_of_nodes = length(mesh.data.NODE);
dofs_per_element = 3;
total_number_of_dofs = total_number_of_nodes * dofs_per_element;

[node_el, etype_element] = mesh.retrievemeshtype(reader);
dofs_per_element = (node_el * dofs_per_element);
sigmaVM0=zeros(total_number_of_nodes,1);

initialdofs = solver.soldofs;
node_elementcorrelation_list = zeros(total_number_of_nodes,2);
for i = 1:total_number_of_elements
    elementTag = mesh_elements(i);
    element_nodes = mesh.data.ELEMENTS{elementTag};
    [nodes,flag,natural_coordinates_all_nodes] = mesh.retrievemeshtype(reader);
    for j = 1:length(element_nodes)
        node_elementcorrelation_list(element_nodes(j),1)=i;
        node_elementcorrelation_list(element_nodes(j),2)=j;
    end
end

parfor i = 1:length(node_elementcorrelation_list(:,1))

    % Recover each element tag
    elementTag = mesh_elements(node_elementcorrelation_list(i,1));
    node_number = node_elementcorrelation_list(i,2);
    element_nodes = mesh.data.ELEMENTS{elementTag};
    number_of_nodes = length(element_nodes);
    element_coordinates=zeros(3,number_of_nodes);
    dof_per_node=2;
    Tee=zeros(1,number_of_nodes);
    Uee=zeros(1,number_of_nodes*3);
    element_dof_indexes=zeros(number_of_nodes*dof_per_node,1);
    element_material_index=mesh.elements_material(elementTag);
    xx= mesh.elements_density(elementTag);
    Number_of_Nodes = length(element_coordinates(1,:));

    for ei=1:number_of_nodes
        element_dof_indexes(ei)=element_nodes(ei)*dof_per_node-1;
        element_dof_indexes(number_of_nodes+i)=element_nodes(ei)*dof_per_node;
        element_coordinates(:,ei)=mesh.data.NODE{element_nodes(ei)};
        Tee(ei)=initialdofs(element_nodes(ei)*dof_per_node-1);
        Uee((ei-1)*3+1)=solver.soldofs_mech((element_nodes(ei)-1)*3+1);
        Uee((ei-1)*3+2)=solver.soldofs_mech((element_nodes(ei)-1)*3+2);
        Uee((ei-1)*3+3)=solver.soldofs_mech((element_nodes(ei)-1)*3+3);
    end
    [nodes,flag,natural_coordinates_all_nodes] = mesh.retrievemeshtype(reader);
            
    j=node_number;
            [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype_element, natural_coordinates_all_nodes(j,1), natural_coordinates_all_nodes(j,2), natural_coordinates_all_nodes(j,3));

            JM = dShape' * element_coordinates';
            %Jacinv = inv(JM);
            DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
            % FIXME, calculate from all dofs input
            Th = N * Tee';

            Dalphapx = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_x');
            Dalphapy = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_y');
            Dalphapz = reader.getmaterialproperty(element_material_index,'ThermalExpansionCoefficient_z');
            DEp = reader.getmaterialproperty(element_material_index,'YoungModulus');
            nu = reader.getmaterialproperty(element_material_index,'PoissonRatio');
            Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_YoungModulus'),reader.getmaterialproperty(element_material_index,'Tmax_YoungModulus')];
            [DE,DdE]=CalculateMaterialProperties(DEp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_YoungModulus'));
            %[Dalpha,Ddalpha]=CalculateMaterialProperties(Dalphap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalExpansionCoefficient'));

            alphav=zeros(6,1);
            alphav(1:3,1)=[Dalphapx,Dalphapy,Dalphapz];

            C = DE./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                0 0 0 0 0 (1-2*nu)/2];

            % Preallocate memory for B-Operators
            Bi=zeros(6,3);B=zeros(6,Number_of_Nodes*3);
            for bi=1:Number_of_Nodes
                Bi(1,1) = DN(1,bi);
                Bi(2,2) = DN(2,bi);
                Bi(3,3) = DN(3,bi);
                Bi(4,:) = [DN(2,bi),DN(1,bi),0];
                Bi(5,:) = [0,DN(3,bi),DN(2,bi)];
                Bi(6,:) = [DN(3,bi),0,DN(1,bi)];
                B(:,(bi-1)*3+1:(bi)*3)=Bi;
            end

            eps=B*Uee'+alphav*N*(Tee-str2double(reader.T0))';
            se0=C*eps-C*alphav*N*(Tee-str2double(reader.T0))';
            sVM0=sqrt(se0(1)^2+se0(2)^2+se0(3)^2-se0(1)*se0(2)-se0(1)*se0(3)-se0(2)*se0(3)+3*se0(4)^2+3*se0(5)^2+3*se0(6)^2);
            sigmaVM0(i) = sVM0;
end
sigmaVM0_ordered=zeros(total_number_of_nodes,1);
for i = 1:length(node_elementcorrelation_list(:,1))
     sigmaVM0_ordered(node_elementcorrelation_list(i,2))=sigmaVM0(i);
end
end

