function [mat,Dmat] = CalculateMaterialProperties(material,Th,xx)
                TC=Th;
                mat=0;Dmat=0;pp=1;
                for i=1:length(material)
                    mat=mat+material(i)*TC^(i-1);
                    if i>1
                        pp=pp*(i-1);
                        Dmat=Dmat+material(i)*TC^(i-2)/pp;
                    end
                end       

end
