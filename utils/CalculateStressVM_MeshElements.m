function [sigmaVM0] = CalculateStressVM_MeshElements(reader,mesh,solver,name)
            meshelements=mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{meshelements(1)};
            dim = mesh.retrieveelementdimension(etype); 
mesh_elements = mesh.retrieveElementalSelection(name);
total_number_of_elements = length(mesh_elements);
total_number_of_nodes = length(mesh.data.NODE);
dofs_per_element = dim;
total_number_of_dofs = total_number_of_nodes * dofs_per_element;

[node_el, etype_element] = mesh.retrievemeshtype(reader);
dofs_per_element = (node_el * dofs_per_element);
sigmaVM0=zeros(total_number_of_elements,1);

initialdofs = solver.soldofs;

for i = 1:total_number_of_elements

    % Recover each element tag
    elementTag = mesh_elements(i);
    element_nodes = mesh.data.ELEMENTS{elementTag};
    number_of_nodes = length(element_nodes);
    element_coordinates=zeros(3,number_of_nodes);
    Tee=zeros(1,number_of_nodes);
    Uee=zeros(1,number_of_nodes*dim);
    element_dof_indexes=zeros(number_of_nodes*2,1);
    element_material_index=mesh.elements_material(elementTag);
    xx= mesh.elements_density(elementTag);
    Number_of_Nodes = length(element_coordinates(1,:));

    for ei=1:number_of_nodes
        element_dof_indexes(ei)=element_nodes(ei)*2-1;
        element_dof_indexes(number_of_nodes+i)=element_nodes(ei)*2;
        element_coordinates(:,ei)=mesh.data.NODE{element_nodes(ei)};
        Tee(ei)=solver.soldofs(element_nodes(ei)*2-1);
        Uee((ei-1)*dim+1)=solver.soldofs_mech((element_nodes(ei)-1)*dim+1);
        Uee((ei-1)*dim+2)=solver.soldofs_mech((element_nodes(ei)-1)*dim+2);
        if dim==3
            Uee((ei-1)*dim+3)=solver.soldofs_mech((element_nodes(ei)-1)*dim+3);
        end
    end
            if dim==2
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, 0,0,0);
                element_coordinates=element_coordinates(1:2,:);
                N=N';
            else
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, 0,0,0);
            end

    %[N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype_element, 0, 0, 0);

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


                if dim==2
               alphav=zeros(3,1);
               alphav(1:2,1)=[Dalphapx,Dalphapy];

                C = DE / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % plane stress
                % C = DE / (1 + nu) / (1 - 2 * nu) * [1 - nu, nu, 0; nu, 1
                % - nu, 0; 0, 0, (1 - 2 * nu) / 2]; % plane strain

                Bi = zeros(3, 2);
                B = zeros(3, Number_of_Nodes * 2);
            
                for bi = 1:Number_of_Nodes
                    Bi(1, 1) = DN(1, bi);
                    Bi(2, 2) = DN(2, bi);
                    Bi(3, 1) = DN(2, bi);
                    Bi(3, 2) = DN(1, bi);
            
                    B(:, (bi - 1) * 2 + 1 : bi * 2) = Bi;
                end                
            elseif dim==3
                alphav=zeros(6,1);
                alphav(1:3,1)=[Dalphapx,Dalphapy,Dalphapz];    

                C = DE./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                    0 0 0 0 0 (1-2*nu)/2];

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
            end

    %eps=B*Uee'+alphav*N*(Tee-str2double(reader.T0))';
    %se0=C*eps-C*alphav*N*(Tee-str2double(reader.T0))';
    %sVM0=sqrt(se0(1)^2+se0(2)^2+se0(3)^2-se0(1)*se0(2)-se0(1)*se0(3)-se0(2)*se0(3)+3*se0(4)^2+3*se0(5)^2+3*se0(6)^2);
    se0=C*B*Uee'-C*alphav*N*(Tee-str2double(reader.T0))';
    %se0=C*eps   -C*alphav*N*(Tee-str2double(reader.T0))';
    if dim==3
        sVM0=sqrt(se0(1)^2+se0(2)^2+se0(3)^2-se0(1)*se0(2)-se0(1)*se0(3)-se0(2)*se0(3)+3*se0(4)^2+3*se0(5)^2+3*se0(6)^2);
    elseif dim==2
        sVM0=sqrt(se0(1)^2+se0(2)^2-se0(1)*se0(2)+3*se0(3)^2);
    end
    sigmaVM0(i) = sVM0;
end

end

