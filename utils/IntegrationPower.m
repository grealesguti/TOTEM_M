function Powe = IntegrationPower(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
                %[N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
            [dim] = mesh.retrieveelementdimension(etype);
            if dim==2
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), 0);
                element_coordinates=element_coordinates(1:2,:);
                N=N';
            else
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));
            end

                JM = dShape' * element_coordinates';
                detJ=det(JM);
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                % FIXME, calculate from all dofs input
                Th = N * Tee';

                Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                Dap = reader.getmaterialproperty(element_material_index,'Seebeck');


                Dep = reader.getmaterialproperty(element_material_index,'ElectricalConductivity');
                Dap = reader.getmaterialproperty(element_material_index,'Seebeck');
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ElectricalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ElectricalConductivity')];
                [De]=CalculateMaterialProperties(Dep,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_Seebeck'),reader.getmaterialproperty(element_material_index,'Tmax_Seebeck')];
                [Da]=CalculateMaterialProperties(Dap,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                %[Dk,Ddk]=CalculateMaterialProperties(Dkp,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                
                je=-De*DN*Vee'-Da*De*DN*Tee';
            
                Powe=detJ*(-je'*DN*Vee');
end

