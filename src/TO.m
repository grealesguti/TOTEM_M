classdef TO
    %TO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m
        n
        xmin
        xmax
        dfdx
        low
        upp
        c
        d
        a0
        a
        outeriter
        maxiter
        inputReader
        mesh
        bcinit
        kkttol
    end
    
    methods
        function obj = TO(inputArg1,inputArg2)
            m=length(con_name);outeriter=0;Conbj=zeros(m,1);
            storagecte=zeros(3+m,1);    storagecte(1:3)=[kkttol Vinf Vsup];
            if dV~=0
            xval=[xx(TOEL)' vx0]';
            n=length(TOEL)+1;
            else
            xval=[xx(TOEL)']';
            n=length(TOEL);
            end
            xmin=zeros(n,1);xmax=ones(n,1);xold1=xval;xold2=xval;
            xmax(end)=scalev;
            dfdx=zeros(m,n);fval=zeros(m,1);
            low     = 0.3;upp     = 0.7;
            c       = ones(m,1)*1000;c(2)=100;
            d       = ones(m,1);
            a0      = 1;
            a       = zeros(m,1);
        end
        
        function outputArg = runMMA(obj,inputArg)
            reader = InputReader("Benchmarks/Elements/input_Benchmark1_LINEARHEX_PARAM.txt");
            fprintf('Initialized InputReader with filename: %s\n', inputfilename);
            mesh = Mesh(reader);
            fprintf('Initialized Mesh\n');
            bcinit = BCInit(reader, mesh);
            solver = Solver();
            
            obj.initMMA() % including dfdx!!! running sensitivities should return dfdx, and filters


        while kktnorm > kkttol && outeriter < maxouteriter 
            outeriter = outeriter+1;
            postprocess.save()
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
            %% Filter densities
            xx=xmma;
            obj.FilteringDensitites()

            %% New NR starting point
            % New Voltage drop
            obj.modifyNRStartingpoint()

            %% New Solve
            solver = Solver();

            %% New derivatives
            dfdx = obj.sensitivities.calculate_dfdx()
            obj.FilteringSensitivities()

            %% MMA parameters update
            xold2=xold1;xold1=xval;
            xval=xmma;

            %% write results
            postprocesing.save()

            %% Convergence
            kktnorm=changexval(outeriter);
            postprocesing.plot()
        end

        end
    end
end

