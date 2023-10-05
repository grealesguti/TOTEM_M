function [mat,Dmat] = CalculateMaterialProperties(material,Th,x,p)
    TC=Th;
    minmat=1e-9;

    mat=0;Dmat=0;pp=1;
    for i=1:length(material)
        mat=mat+material(i)*TC^(i-1);
        if i>1
            pp=pp*(i-1);
            Dmat=Dmat+material(i)*TC^(i-2)*(i-1);
        end
    end       

    mat=minmat+x^p*(mat-minmat);
    Dmat=Dmat*x^p;

end
