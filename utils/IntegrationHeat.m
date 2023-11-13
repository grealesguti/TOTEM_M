function Qth = IntegrationHeat(natural_coordinates, element_coordinates, Tee, Vee, element_material_index, reader, mesh, etype,xx)
                [N, dShape] = mesh.selectShapeFunctionsAndDerivatives(etype, natural_coordinates(1), natural_coordinates(2), natural_coordinates(3));

                JM = dShape' * element_coordinates';
                detJ=det(JM);
                %Jacinv = inv(JM);
                DN = inv(JM) * dShape'; % FIXME and check it is the same! NOT THE SAME RESULT!!!
                % FIXME, calculate from all dofs input
                Th = N * Tee';

                Dkp = reader.getmaterialproperty(element_material_index,'ThermalConductivity');
                %[De]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_ElectricalConductivity'));
                %[Da]=CalculateMaterialProperties(Dap,Th,xx,reader.getmaterialproperty(element_material_index,'Penalty_Seebeck'));
                Tmat=[Th,reader.getmaterialproperty(element_material_index,'Tmin_ThermalConductivity'),reader.getmaterialproperty(element_material_index,'Tmax_ThermalConductivity')];
                [Dk,Ddk]=CalculateMaterialProperties(Dkp,Tmat,xx,reader.getmaterialproperty(element_material_index,'Penalty_ThermalConductivity'));
                Qth=detJ*(-Dk*DN*Tee);
end

