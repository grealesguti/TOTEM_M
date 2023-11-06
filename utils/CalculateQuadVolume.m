function[vt]=CalculateQuadVolume(nodeposh)
        % numerically calculates the volume of any hexahedral given the
        % nodal coordinates of its 8 vertex 
        
        % nodal  positions of each triangle
        th(1,:)=[1, 2, 3, 5];
        th(2,:)=[4, 2, 6, 5];
        th(3,:)=[4, 6, 8, 5];
        th(4,:)=[4, 2, 3, 7];
        th(5,:)=[4, 6, 2, 7];
        th(6,:)=[4, 8, 6, 7];
        % nodar order of each tetrahedra face
        faces2(1,:)=[1, 3, 2];
        faces2(2,:)=[2, 4, 1];
        faces2(3,:)=[3, 4, 2];
        faces2(4,:)=[3, 1, 4];
        
        vt=0.;      % initialization of total volume
        for k=1:6   % loop over tetrahedra
              nodepos(1:4,:)=nodeposh(th(k,1:4),:)  ; %
              V=0.;
            for i =1:4     % loop over tetrahedra faces
                    fvertex(1,:)=nodepos(faces2(i,1),:);
                    fvertex(2,:)=nodepos(faces2(i,2),:);
                    fvertex(3,:)=nodepos(faces2(i,3),:);

                    xi=fvertex(:,1);
                    yi=fvertex(:,2);
                    zi=fvertex(:,3);
                    
                beta=[yi(2)-yi(1),yi(3)-yi(1)];
                gam=[zi(2)-zi(1),zi(3)-zi(1)];
                d0=beta(1)*gam(2)-beta(2)*gam(1);

                V=V+d0*sum(xi)/6.;
                
            end

            cmv(k,1)=sum(nodeposh(:,1))/4.;
            cmv(k,2)=sum(nodeposh(:,2))/4.;
            cmv(k,3)=sum(nodeposh(:,3))/4.;
            vv(k)=V;
            vt=V+vt;
        end

        %cmf(1)=(cmv(:,1)*vv)/vt/2.
        %cmf(2)=(cmv(:,2)*vv)/vt/2.
        %cmf(3)=(cmv(:,3)*vv)/vt/2.
        cmf=[];
        end