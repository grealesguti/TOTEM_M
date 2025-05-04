function [mat_x,Dmat_x] = CalculateMaterial_XDerivative(minmat,material,T,x,p)
    TC=T(1);
    Tmin=T(2);
    Tmax=T(3);

    if TC<Tmin
        TC=Tmin;
    elseif TC>Tmax
        TC=Tmax;
    end    

    %minmat=1e-6;
%     mod=0;
%     if TC<325
%         TC=325;
%         mod=1;
%     elseif TC>600
%         TC=600;
%         mod=1;
%     end
    mat=0;Dmat=0;pp=1;
    for i=1:length(material)
        mat=mat+material(i)*TC^(i-1);
        if i>1
            pp=pp*(i-1);
            Dmat=Dmat+material(i)*TC^(i-2)*(i-1);
        end
    end       

    mat_x=p*x^(p-1)*(mat-minmat);
    Dmat_x=Dmat*p*x^(p-1);
    %Dmat_x=0; % parcial derivative, no degrees of freedom, including T!!!
    if TC<Tmin
        Dmat_x=0;
    elseif TC>Tmax
        Dmat_x=0;
    end
end
