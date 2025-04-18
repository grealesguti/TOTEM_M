classdef Postprocessing < handle
    %POSTPROCESSING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        total_number_of_elements
        total_number_of_nodes
        coordinates
        element_node_idxs
        mesh_elements
        dim
    end
    
    methods
        function obj = Postprocessing()
       
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initVTK(obj,reader,mesh)
            obj.mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            obj.total_number_of_elements = length(obj.mesh_elements);
            obj.total_number_of_nodes = length(mesh.data.NODE);
                
            obj.coordinates = zeros(3,obj.total_number_of_nodes);
    
            for i=1:obj.total_number_of_nodes
                obj.coordinates(:,i)=mesh.data.NODE{i};
            end
            meshelements=mesh.retrieveElementalSelection(reader.MeshEntityName);
            etype=mesh.data.ElementTypes{meshelements(1)};
            obj.dim = mesh.retrieveelementdimension(etype); 
            if obj.dim==2
                node_per_el = 4;
            else
                node_per_el = 8;
            end

            obj.element_node_idxs=zeros(obj.total_number_of_elements,node_per_el);
            for i=1:obj.total_number_of_elements
                nodal_idxx=mesh.data.ELEMENTS{obj.mesh_elements(i)};
                obj.element_node_idxs(i,:) = nodal_idxx(1:node_per_el);
            end        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_Mesh(obj,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_Mesh_', dateStr, '.vtk');         
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs,[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_Mesh_xx(obj,filepath,mesh)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_Mesh_', dateStr, '.vtk');         
       
            xx = mesh.elements_density;

            printxx=xx(obj.mesh_elements)' ;
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'CELL_DATA',printxx,...
                'data',[],'precision',5)      

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_U(obj,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            %dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            %outputFilePath = append(filepath, '_U_', dateStr, '.vtk');         
            outputFilePath=filepath;
            Ux=zeros(obj.total_number_of_nodes,1);
            Uy=zeros(obj.total_number_of_nodes,1);
            Uz=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Ux(i) = solver.soldofs_mech((i-1)*obj.dim+1);
                Uy(i) = solver.soldofs_mech((i-1)*obj.dim+2);
                if obj.dim==3
                    Uz(i) = solver.soldofs_mech((i-1)*obj.dim+obj.dim);
                end
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',[Ux';Uy';Uz']','U','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_U_axis(obj,solver,axis,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_U_', dateStr, '.vtk');         

            Ux=zeros(obj.total_number_of_nodes,1);
            Uy=zeros(obj.total_number_of_nodes,1);
            Uz=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Ux(i) = solver.soldofs_mech((i-1)*3+1);
                Uy(i) = solver.soldofs_mech((i-1)*3+2);
                Uz(i) = solver.soldofs_mech((i-1)*3+3);
            end
            if axis==1
                U=Ux;
            elseif axis==2
                U=Uy;
            elseif axis==3
                U=Uz;
            end
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',U,append(['U',num2str(axis)]),'Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_TV(obj,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');         

            Tn=zeros(obj.total_number_of_nodes,1);
            Vn=zeros(obj.total_number_of_nodes,1);
            
            for i=1:obj.total_number_of_nodes
                Tn(i) = solver.soldofs(i*2-1);
                Vn(i) = solver.soldofs(i*2);
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',[Tn';Vn']','TV','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_qth(obj,reader,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');         

            qth = CalculateHeatFlow_nodes(reader,mesh,solver);
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',qth','qth','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_q(obj,reader,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');         

            q = CalculateTotalHeatFlow_nodes(reader,mesh,solver);
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',q','qth','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_j(obj,reader,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');         

            j = CalculateTotalCurrentDensity_nodes(reader,mesh,solver);
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'data','POINT_DATA',j','qth','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_x_TV(obj,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            outputFilePath = filepath;         
            
            xx = mesh.elements_density;

            Tn=zeros(obj.total_number_of_nodes,1);
            Vn=zeros(obj.total_number_of_nodes,1);
            
            for i=1:obj.total_number_of_nodes
                Tn(i) = solver.soldofs(i*2-1);
                Vn(i) = solver.soldofs(i*2);
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'CELL_DATA',xx(obj.mesh_elements)',...
                'data','POINT_DATA',[Tn';Vn']','TV','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_x_U(obj,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            outputFilePath = filepath;         
            
            xx = mesh.elements_density;

            % Append the date and '.vtk' extension to the filepath
            %dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            %outputFilePath = append(filepath, '_U_', dateStr, '.vtk');         
            
            Ux=zeros(obj.total_number_of_nodes,1);
            Uy=zeros(obj.total_number_of_nodes,1);
            Uz=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Ux(i) = solver.soldofs_mech((i-1)*obj.dim+1);
                Uy(i) = solver.soldofs_mech((i-1)*obj.dim+2);
                if obj.dim==3
                    Uz(i) = solver.soldofs_mech((i-1)*obj.dim+obj.dim);
                end
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'CELL_DATA',xx(obj.mesh_elements)',...
                'data','POINT_DATA',[Ux';Uy';Uz']','U','Test',[],'precision',5)             
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_matidx_TV(obj,mesh,solver,filepath)
            % Append the date and '.vtk' extension to the filepath
            %outputFilePath = filepath;         
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');     
            xx = mesh.elements_material;

            Tn=zeros(obj.total_number_of_nodes,1);
            Vn=zeros(obj.total_number_of_nodes,1);
            
            for i=1:obj.total_number_of_nodes
                Tn(i) = solver.soldofs(i*2-1);
                Vn(i) = solver.soldofs(i*2);
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'CELL_DATA',xx(obj.mesh_elements)',...
                'material','POINT_DATA',[Tn';Vn']','TV','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function VTK_freedofs(obj,bcinit,filepath)
            % Append the date and '.vtk' extension to the filepath
            %outputFilePath = filepath;         
            dateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
            outputFilePath = append(filepath, '_TV_', dateStr, '.vtk');     
            fixdofs = bcinit.dofs_fixed_;

            Tn=zeros(obj.total_number_of_nodes,1);
            Vn=zeros(obj.total_number_of_nodes,1);
            
            for i=1:obj.total_number_of_nodes
                Tn(i) = fixdofs(i*2-1);
                Vn(i) = fixdofs(i*2);
            end
    
            vtkwrite( outputFilePath, ...
                'unstructured_grid',obj.coordinates(1,:),obj.coordinates(2,:),obj.coordinates(3,:),...
                'CELLS',obj.element_node_idxs, ...
                'POINT_DATA',[Tn';Vn']','TVfixed','Test',[],'precision',5)        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function PlotIter(~,fig,reader,iter,f0val,fval,xbc)
            % Ensure the specified figure exists and is active or create a new one
            if isempty(fig) || ~ishandle(fig) || ~strcmp(get(fig, 'Type'), 'figure')
                fig = figure;
                set(fig, 'Position', [100, 100, 1200, 400]);
            else
                figure(fig);
            end
        
            subplot(1, 1 + length(fval(1, :))+length(reader.TObcval),1);
            % Plot the objective function history
            plot(1:iter, f0val(1:iter), 'o-');
            title(reader.TopOpt_Objective);
        
            % Plot constraint histories in subplots
            for i = 1:length(fval(1, :))
                subplot(1, 1 + length(fval(1, :))+length(reader.TObcval),i+1);
                hold on
                yyaxis left;
                plot(1:iter, fval(1:iter, i), 'o-');
                yline(0,'b')
                title(reader.TopOpt_ConstraintName{i});
                yyaxis right;
                plot(1:iter, (fval(1:iter, i)+1)*reader.TopOpt_ConstraintValue(i), 'g--');
                yline(reader.TopOpt_ConstraintValue(i),'r--')

            end

            if not(isempty(reader.TObcval))
                for i = 1:length(reader.TObcval)
                    subplot(1, 1 + length(fval(1, :))+length(reader.TObcval),i+1+length(fval(1, :)));
                    hold on
                    xbc_value=xbc(1:iter,i);
                    bc_value= reader.TObcminval(i)+xbc_value*(reader.TObcmaxval(i)-reader.TObcminval(i));
                    yyaxis left;
                    plot(1:iter, xbc_value, 'o-');
                    yline(1,'b')
                    yline(0,'b')
                    title(reader.TObctype{i});
                    yyaxis right;
                    plot(1:iter, bc_value, 'g--');
                    yline(reader.TObcmaxval(i),'r--')
                    yline(reader.TObcminval(i),'r--')
                end
            end
            %folder = reader.rst_folder;
                        % Split the string by '/' symbol
            %splitString = strsplit(folder, '/');
            
            % Recover the last element
            %lastElement = splitString{end};
            % Update the figure
            sgtitle(append('Optimization Progress:',reader.Rst_name));
            hold off;

        % Add this line to update the figure in each iteration
        drawnow;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveIterCSV(~, filepath, reader, iter, f0val, fval, xbc, kktnorm)
    
    % Initialize a cell array to store the data
    dataToWrite = cell(iter, 3 + length(fval(1, :)) + length(reader.TObcval));
    
    % Add headers to the cell array
    dataToWrite{1, 1} = 'Iteration';
    dataToWrite{1, 2} = reader.TopOpt_Objective;
    constraintNames = reader.TopOpt_ConstraintName;
    
    for i = 1:length(fval(1, :))
        dataToWrite{1, 2 + 2 * i - 1} = [constraintNames{i}, '_value'];
        dataToWrite{1, 2 + 2 * i} = [constraintNames{i}, '_right_value'];
    end
    
    if ~isempty(reader.TObcval)
        for i = 1:length(reader.TObcval)
            dataToWrite{1, 2 + 2 * length(fval(1, :)) + i} = [reader.TObctype{i}, '_value'];
            dataToWrite{1, 2 + 2 * length(fval(1, :)) + i + 1} = [reader.TObctype{i}, '_right_value'];
        end
    end
    
    % Add header for kktnorm
    dataToWrite{1, end} = 'KKTNorm';
    
    % Store data for each iteration
    for i = 2:iter
        dataToWrite{i, 1} = i - 1;
        dataToWrite{i, 2} = f0val(i);
        for j = 1:length(fval(1, :))
            dataToWrite{i, 2 + 2 * j - 1} = fval(i, j);
            dataToWrite{i, 2 + 2 * j} = (fval(i, j) + 1) * reader.TopOpt_ConstraintValue(j);
        end
        if ~isempty(reader.TObcval)
            for k = 1:length(reader.TObcval)
                dataToWrite{i, 2 + 2 * length(fval(1, :)) + k} = xbc(i, k);
                dataToWrite{i, 2 + 2 * length(fval(1, :)) + k + 1} = reader.TObcminval(k) + xbc(i, k) * (reader.TObcmaxval(k) - reader.TObcminval(k));
            end
        end
        % Store the kktnorm value for this iteration
        dataToWrite{i, end} = sprintf('%.3e', kktnorm(i)); % 3 decimal places in scientific notation
    end
    
    % Specify the CSV file path
    csvFilePath = filepath;
    
    % Write the data to a CSV file
    writecell(dataToWrite, csvFilePath, 'WriteMode', 'append');
end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sorted_Tn]=Benchmark_T_PLOT_axis(obj,fig,solver,axis)
            Tn=zeros(obj.total_number_of_nodes,1);
            Tn_loc=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Tn(i) = solver.soldofs(i*2-1);
                Tn_loc(i) = obj.coordinates(axis,i);
            end
            % Sort Tn_loc and get the sorting indices
            [sorted_locs, sorting_indices] = sort(Tn_loc);
            
            % Rearrange Tn and Tn_loc based on the sorting indices
            sorted_Tn = Tn(sorting_indices);
            sorted_Tn_loc = Tn_loc(sorting_indices);
            
            %figure(fig)
            % Set the background color of the figure to white
            set(gcf, 'Color', 'white')
            plot(sorted_Tn_loc, sorted_Tn, '-o')  % Use 'o' for markers
            ylabel('Temperature [K]')
            xlabel('Location [m]')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sorted_Vn]=Benchmark_V_PLOT_axis(obj,fig,solver,axis)

            Vn=zeros(obj.total_number_of_nodes,1);
            Vn_loc=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Vn(i) = solver.soldofs(i*2);
                Vn_loc(i) = obj.coordinates(axis,i);
            end
            % Sort Tn_loc and get the sorting indices
            [sorted_locs, sorting_indices] = sort(Vn_loc);
            
            % Rearrange Tn and Tn_loc based on the sorting indices
            sorted_Vn = Vn(sorting_indices);
            sorted_Vn_loc = Vn_loc(sorting_indices);
            %figure(fig)
            plot(sorted_Vn_loc,sorted_Vn, '-o')
            ylabel('Voltage [V]')
            xlabel('Location [m]')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sorted_Un]=Benchmark_U_PLOT_axis(obj,fig,solver,axis)
            Un=zeros(obj.total_number_of_nodes,1);
            Ux=zeros(obj.total_number_of_nodes,1);
            Uy=zeros(obj.total_number_of_nodes,1);
            Uz=zeros(obj.total_number_of_nodes,1);

            Un_loc=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Ux(i)= solver.soldofs_mech((i-1)*obj.dim+1);
                Uy(i)= solver.soldofs_mech((i-1)*obj.dim+2);
                if obj.dim==3
                    Uz(i)= solver.soldofs_mech((i-1)*3+3);
                end
                Un(i) = solver.soldofs_mech((i-1)*obj.dim+axis);
                Un_loc(i) = obj.coordinates(axis,i);
            end
            % Sort Tn_loc and get the sorting indices
            [sorted_locs, sorting_indices] = sort(Un_loc);
            
            % Rearrange Tn and Tn_loc based on the sorting indices
            sorted_Un = Un(sorting_indices);
            sorted_Un_loc = Un_loc(sorting_indices);
            
            %figure(fig)
            % Set the background color of the figure to white
            set(gcf, 'Color', 'white')
            plot(sorted_Un_loc, sorted_Un, '-o')  % Use 'o' for markers
            ylabel('Displacement [m]')
            xlabel('Location [m]')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotmaterialT(obj,fig,element_material_index,material_property,Tmin,Tmax,reader)
            steps=100;
            Dev=zeros(1,steps);
            Ddev=zeros(1,steps);
            Thv=zeros(1,steps);
            c=1;
            xx=1;
            for Th = Tmin:(Tmax-Tmin)/steps:Tmax
                Dep = reader.getmaterialproperty(element_material_index,material_property);
                [De,Dde]=CalculateMaterialProperties(Dep,Th,xx,reader.getmaterialproperty(element_material_index,append('Penalty_',material_property)));
                Dev(c)=De;
                Ddev(c)=Dde;
                Thv(c)=Th;
                c=c+1;
            end

            figure(fig)
            yyaxis left
            plot(Thv,Dev)
            yyaxis right
            plot(Thv,Ddev)
            grid on
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end
end

