classdef Postprocessing < handle
    %POSTPROCESSING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        total_number_of_elements
        total_number_of_nodes
        coordinates
        element_node_idxs
        mesh_elements
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
            obj.element_node_idxs=zeros(obj.total_number_of_elements,8);
            for i=1:obj.total_number_of_elements
                nodal_idxx=mesh.data.ELEMENTS{obj.mesh_elements(i)};
                obj.element_node_idxs(i,:) = nodal_idxx(1:8);
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
        function VTK_U(obj,solver,filepath)
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
        
            % Update the figure
            sgtitle('Optimization Progress');
            hold off;

        % Add this line to update the figure in each iteration
        drawnow;
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
                Ux(i)= solver.soldofs_mech((i-1)*3+1);
                Uy(i)= solver.soldofs_mech((i-1)*3+2);
                Uz(i)= solver.soldofs_mech((i-1)*3+3);
                Un(i) = solver.soldofs_mech((i-1)*3+axis);
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
            ylabel('Temperature [K]')
            xlabel('Location [m]')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end
end

