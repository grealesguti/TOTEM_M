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
             aa=mesh.retrieveElementalSelection(reader.MeshEntityName);
            if (all==1)
                mesh.elements_density(aa)=xx';
            elseif (all<1)
                for i = 1: length(xx)
                    if xx(i)>all
                        mesh.elements_density(aa(i))=1;
                    else
                        mesh.elements_density(aa(i))=1e-5;
                    end
                end
            else
                TOEL=mesh.retrieveElementalSelection(reader.TopOpt_DesignElements);

                mesh.elements_density(TOEL)=ones(length(mesh.elements_density(TOEL)),1)*xx;
            end

end

