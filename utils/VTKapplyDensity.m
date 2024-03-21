function [] = VTKapplyDensity(reader,mesh,all)
            filename = reader.TopOpt_Initial_x;

            A=regexp(fileread((filename)),'\n','split');
            An=find(startsWith(A,'CELL_DATA')==1);
            B=strtrim(strsplit(char(A(startsWith(A,'CELL_DATA')==1)),' '));
            Nxxdata=str2double(char(B(2)));
            xx=zeros(Nxxdata,1);
            for i=1:Nxxdata
                xx(i)=str2double(char(A(An+2+i)));
            end

            if (all==1)
                    mesh.elements_density=xx';
            else
                TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);

                mesh.elements_density(TOEL)=ones(length(mesh.elements_density(TOEL)),1)*xx;
            end

end

