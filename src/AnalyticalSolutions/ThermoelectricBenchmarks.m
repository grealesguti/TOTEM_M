classdef ThermoelectricBenchmarks
    %THERMOELECTRICBENCHMARKS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Benchmarks
        BenchmarksFunctions
    end
    
    methods
        function obj = ThermoelectricBenchmarks()
            close all;
            %THERMOELECTRICBENCHMARKS Construct an instance of this class
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear = "Benchmarks/Elements/Benchmark_HexLinear/input_LinearUncoupledSeebeck.txt";
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear = "Benchmarks/Elements/Benchmark_HexLinear/input_NonLinCouplSEffect.txt";
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity = "Benchmarks/Elements/Benchmark_HexSerendipity/input_LinearUncoupledSeebeck.txt";
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity = "Benchmarks/Elements/Benchmark_HexSerendipity/input_NonLinCouplSEffect.txt";

            obj.Benchmarks={
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexLinear;
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexLinear;
            Benchmark_Perez_Aparicio_LinUncoupSEffect_HexSerendipity;
            Benchmark_Perez_Aparicio_NonLinCouplSEffect_HexSerendipity};

            obj.BenchmarksFunctions=[1,2,1,2];

            for i=1:length(obj.BenchmarksFunctions)
                if(obj.BenchmarksFunctions(i)==1)
                    Perez_Aparicio_LinearUncoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==2)
                    Perez_Aparicio_CoupledSeebeck(obj,obj.Benchmarks{i},i)
                elseif(obj.BenchmarksFunctions(i)==3)
                    Perez_Aparicio_CoupledPeltier(obj,obj.Benchmarks{i},i)
                end

            end
        end
        % https://link.springer.com/article/10.1007/s00466-006-0080-7
        function [xv,Tx,Vx] = Perez_Aparicio_LinearUncoupledSeebeck(~,filename,index)

                reader = InputReader(filename);
                fprintf('Initialized InputReader with filename: %s\n', filename);
                mesh = Mesh(reader);
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                postprocessing = Postprocessing();
                postprocessing.initVTK(reader,mesh);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                figure(index)
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=412-sqrt(169587-1.31e7*xv(i))+273.15;
                    Vx(i)=sqrt(0.018-1.31*xv(i))+2.36*xv(i)-0.13;
                end
                plot(xv,Tx)

        end

        function [xv,Tx,Vx] = Perez_Aparicio_CoupledSeebeck(~,filename,index)

                reader = InputReader(filename);
                fprintf('Initialized InputReader with filename: %s\n', filename);
                mesh = Mesh(reader);
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                postprocessing = Postprocessing();
                postprocessing.initVTK(reader,mesh);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                figure(index)
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=3.794e7*xv(i)*(1.524e-3-xv(i))+298.15;
                    Vx(i)=5.788e-2-4.913*10*xv(i)+7.315e3*xv(i)^2;
                end

                plot(xv,Tx)

        end        

        function [xv,Tx,Vx] = Perez_Aparicio_CoupledPeltier(~,filename,index)

                reader = InputReader(filename);
                fprintf('Initialized InputReader with filename: %s\n', filename);
                mesh = Mesh(reader);
                fprintf('Initialized Mesh\n');
                bcinit = BCInit(reader, mesh);
                fprintf('Initialized Loads\n');
                solver = Solver(mesh, bcinit);
                solver.runNewtonRaphson(reader, mesh, bcinit);
                fprintf('Postprocessing\n');
                postprocessing = Postprocessing();
                postprocessing.initVTK(reader,mesh);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                figure(index)
                hold on
                postprocessing.Benchmark_T_PLOT_axis(index,solver,2)
                L=1.524e-3;
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tx=zeros(steps,1);
                Vx=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tx(i)=3.794e7*xv(i)*(1.524e-3-xv(i))+298.15;
                    Vx(i)=5.788e-2-4.913*10*xv(i)+7.315e3*xv(i)^2;
                end

                plot(xv,Tv)

        end        
        function [xv,Tv,Vv] = Benchmark1_2_3(~,reader,mesh,solver,postprocessing)
            %METHOD1 Summary of this method goes here
                g1 = reader.getmaterialproperty(1,'ElectricalConductivity');
                k1 = reader.getmaterialproperty(1,'ThermalConductivity');
                a1 = reader.getmaterialproperty(1,'Seebeck');    
                coordinates=zeros(3,length(mesh.data.NODE));
                for i=1:length(mesh.data.NODE)
                    coordinates(:,i)=mesh.data.NODE{i};
                end
                L=abs(max(coordinates(2,:))-min(coordinates(2,:)));
                l1=abs(max(coordinates(1,:))-min(coordinates(1,:)));
                l2=abs(max(coordinates(3,:))-min(coordinates(3,:)));
                A=l1*l2;
                boundaryConditions = reader.bctype;
                % Iterate through the cell array to find the index location of the string
                for i = 1:numel(boundaryConditions)
                    if strcmp(boundaryConditions{i}, 'heat_n')
                        qin = reader.bcval(i);
                    elseif strcmp(boundaryConditions{i}, 'Temperature')
                        Th = reader.bcval(i);
                    elseif strcmp(boundaryConditions{i}, 'Voltage')
                        if(abs(reader.bcval(i))>0)
                            Vc=reader.bcval(i);
                        end
                    end
                end
                V0=0;
                syms C31 Tc j1 
                eqTc1 = -j1^2 / (k1*g1) * L^2/2 + C31*L + Th == Tc;
                eqVc1 = -j1/g1*L + a1*j1^2/(g1*k1)*L^2/2 - a1*C31*L == Vc;              
                eqq = (a1*j1*Tc-k1*(-j1^2/(k1*g1)*L+C31))*A == -qin;
            
                Ss=[eqTc1, eqVc1, eqq];
                % SoLve the equations
                S = solve(Ss, ...
                         [Tc, C31, j1]); 
                TCs=eval(S.Tc)
                C31s=eval(S.C31)
                Vcs=eval(S.Vc)
                steps=100;
                x_step=L/steps;
                xv=zeros(steps,1);
                Tv=zeros(steps,1);
                Vv=zeros(steps,1);
                for i=1:100
                    xv(i)=x_step*(i-1);
                    Tv(i)=-j1^2 / (k1*g1) * xv(i)^2/2 + C31*xv(i) + Th;
                    Vv(i)=-j1/g1*xv(i) + a1*j1^2/(g1*k1)*xv(i)^2/2 - a1*C31*xv(i) + V0;
                end

                figure(1)
                plot(xv,Tv)

        end
    end
end

