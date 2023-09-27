classdef Postprocessing < handle
    %POSTPROCESSING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        total_number_of_elements
        total_number_of_nodes
        coordinates
        element_node_idxs
    end
    
    methods
        function obj = Postprocessing()
       
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initVTK(obj,reader,mesh)
            mesh_elements = mesh.retrieveElementalSelection(reader.MeshEntityName);
            obj.total_number_of_elements = length(mesh_elements);
            obj.total_number_of_nodes = length(mesh.data.NODE);
                
            obj.coordinates = zeros(3,obj.total_number_of_nodes);
    
            for i=1:obj.total_number_of_nodes
                obj.coordinates(:,i)=mesh.data.NODE{i};
            end
            obj.element_node_idxs=zeros(obj.total_number_of_elements,8);
            for i=1:obj.total_number_of_elements
                nodal_idxx=mesh.data.ELEMENTS{mesh_elements(i)};
                obj.element_node_idxs(i,:) = nodal_idxx(1:8);
            end        
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
        function Benchmark_T_PLOT_axis(obj,fig,solver,axis)
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
            
            figure(fig)
            % Set the background color of the figure to white
            set(gcf, 'Color', 'white')
            plot(sorted_Tn_loc, sorted_Tn, '-o')  % Use 'o' for markers
            ylabel('Temperature [K]')
            xlabel('Location [m]')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Benchmark_V_PLOT_axis(obj,fig,solver,axis)
            Vn=zeros(obj.total_number_of_nodes,1);
            Vn_loc=zeros(obj.total_number_of_nodes,1);
            for i=1:obj.total_number_of_nodes
                Vn(i) = solver.soldofs(i*2-1);
                Vn_loc(i) = obj.coordinates(axis,i);
            end
            
            figure(fig)
            plot(Vn_loc,Vn)
            ylabel('Voltage [V]')
            xlabel('Location [m]')
        end
    end
end

