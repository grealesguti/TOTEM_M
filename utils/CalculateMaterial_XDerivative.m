function [mat_x,Dmat_x] = CalculateMaterial_XDerivative(material,Th,x,p)
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

    mat_x=x^(p-1)*(mat-minmat);
    Dmat_x=Dmat*x^(p-1);
    Dmat_x=0;

end
