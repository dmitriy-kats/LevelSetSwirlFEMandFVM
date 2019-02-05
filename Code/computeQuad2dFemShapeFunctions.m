function[xqR, wqR, NR, dNdxR] = computeQuad2dFemShapeFunctions(coords)

x=coords(:,1); %x coordinates
y=coords(:,2); %y coordinates

    %quad points
    zeta=1/sqrt(3).*[-1 1 1 -1]; 
    eta=1/sqrt(3).*[-1 -1 1 1];

    for q=1:4 % loop over quad points
    N1=1/4*(1-zeta(q))*(1-eta(q));
    N2=1/4*(1+zeta(q))*(1-eta(q));
    N3=1/4*(1+zeta(q))*(1+eta(q));
    N4=1/4*(1-zeta(q))*(1+eta(q));
    N(q,:)=[N1 N2 N3 N4]; %shape functions
    %shape funciton derivatives with zeta/eta
    dNdzeta(q,:)=1/4.*[-(1-eta(q)) (1-eta(q)) (1+eta(q)) -(1+eta(q))];
    dNdeta(q,:)=1/4.*[-(1-zeta(q)) -(1+zeta(q)) (1+zeta(q)) (1-zeta(q))];

    xq(q)=[N1 N2 N3 N4]*x;
    yq(q)=[N1 N2 N3 N4]*y;

     dxdzeta(q)=dNdzeta(q,:)*x;
     dxdeta(q)=dNdeta(q,:)*x;
     dydzeta(q)=dNdzeta(q,:)*y;
     dydeta(q)=dNdeta(q,:)*y;
     Jlong(q)=dxdzeta(q)*dydeta(q)-dydzeta(q)*dxdeta(q); %weights
     
     for p=1:4 %shape funciton derivatives with x/y
     tempN=1/(dxdzeta(q)*dydeta(q)-dydzeta(q)*dxdeta(q))...
         *([ dydeta(q) -dydzeta(q); -dxdeta(q) dxdzeta(q)])*[dNdzeta(q,p) dNdeta(q,p)]';
     dNdx(q,p)=tempN(1);
     dNdy(q,p)=tempN(2);
     end      
        
    end
    
    %rearrange to make it work with the rest of the code
    J=reshape(Jlong, 2, 2);   
    xqR=[xq' yq'];
    wqR=Jlong';
    NR=N;
    dNdxR(:,:,1)=dNdx;
    dNdxR(:,:,2)=dNdy;
    
    

end



