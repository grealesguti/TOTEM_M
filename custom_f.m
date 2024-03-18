function loss = custom_f()

        filepath="TECTO/input_TECTO_Thermoel_Serend_cte.txt";   
        reader = InputReader(filepath);
        reader.Rst_name=append('Al2O3_nonlin_0topU_3D_Q_',num2str(qinval),'_P',num2str(powval));

        mesh = Mesh(reader);
        TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);
            bcinit = BCInit(reader, mesh);
                if reader.Filter>0
                    filtering = Filtering(reader,mesh);
                end

                %% Filter densities
                for i = 1:length(reader.TObcval)
                    if(strcmp(reader.TObctype,'Voltage'))
                        Voltage_value=reader.TObcminval(i)+x(length(TOEL)+i)*(reader.TObcmaxval(i)-reader.TObcminval(i));
                        reader.TObcval(i)=Voltage_value;
                        reader.bcval(length(reader.bcval)-length(reader.TObcval)+i)=Voltage_value;
                    end
                end
                %mesh_1 = Mesh(reader);
                for i=1:length(obj.TOEL)
                    mesh.elements_density(obj.TOEL(i))=x(i);
                end
                if reader.Filter>0
                    filtering.filter_densities(reader,mesh)
                end

                %bcinit1 = BCInit(reader, mesh);
                solver = Solver(mesh, bcinit);

                for i=1:length(reader.TObcval)
                    nodes=mesh.retrieveNodalSelection(reader.TObcloc(i));
                    if(strcmp(reader.TObctype,'Voltage'))
                        %solver.soldofs(odd_numbers)=prevdofs_odd;
                        %solver.soldofs(even_numbers)=prevdofs_even*obj.Voltage_value/obj.Voltage_initial;
                        solver.soldofs(nodes*2)=obj.Voltage_value;
                    end
                end

            if strcmp(reader.solver,'NR')
                    residual_norm=solver.runNewtonRaphson(reader, mesh, bcinit);
                    odd_numbers = 1:2:length(solver.soldofs);
                    Tdofs=solver.soldofs(odd_numbers);
                    if residual_norm>10000  || not(isempty(Tdofs(Tdofs<0)))% divergence in NR catch
                        warning('NR DIVERGED, changing to Arc-len!!!');
                        for i=1:length(bcinit.dofs_free_)
                            df=bcinit.dofs_free_(i);
                            if mod(df, 2)==0
                                solver.soldofs(df)=0.0;
                            else 
                                solver.soldofs(df)=0;
                            end
                        end
                        funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                        [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                        solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                        if strcmp(reader.physics,'decoupledthermoelectromechanical')
                            % Extract the necessary variables
                            [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                            Temperature_solution = solver.soldofs(1:2:end);
                            solver.KUT=KThermalLoad;
                            solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                            solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                        end
                    end
                elseif strcmp(reader.solver,'Arc-len')
                        for i=1:length(solver.soldofs)/2
                            if solver.soldofs(i*2-1)==0
                                solver.soldofs(i*2-1)=str2double(reader.T0);
                            elseif solver.soldofs(i*2)==0
                                solver.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                else
                    warning('No solver recognized, changing to Arc-len!!!');
                        for i=1:length(solver.soldofs)/2
                            if solver.soldofs(i*2-1)==0
                                solver.soldofs(i*2-1)=str2double(reader.T0);
                            end
                            if solver.soldofs(i*2)==0
                                solver.soldofs(i*2)=0.01;
                            end
                        end
                    funAL = @(t) solver.funArcLen(reader,bcinit,mesh,t);
                    [ufree] = arc_length_Crisfield(funAL,solver.soldofs(bcinit.dofs_free_));
                    solver.soldofs(bcinit.dofs_free_)=ufree(:,end);
                    if strcmp(reader.physics,'decoupledthermoelectromechanical')
                        % Extract the necessary variables
                        [solver.KStiff,KThermalLoad] = solver.Assembly_DecoupledThermoMech(reader, mesh);
                        Temperature_solution = solver.soldofs(1:2:end);
                        solver.KUT=KThermalLoad;
                        solver.loadVector_mech=solver.loadVector_mech+KThermalLoad*(Temperature_solution-str2double(reader.T0));
                        solver.SolveLinearSystemInParallel(reader.physics,bcinit)
                    end
                end

            TOO = TO_Objectives(reader,mesh,bcinit);
            TOO.CalculateObjective(reader,mesh,solver)
            TOC = TO_Constraints(reader,mesh,bcinit);
            TOC.CalculateConstraint(reader,mesh,solver);

            %% New derivatives
            f0val = TOO.fval;
            fval = TOC.fval;

            loss=f0val;
end