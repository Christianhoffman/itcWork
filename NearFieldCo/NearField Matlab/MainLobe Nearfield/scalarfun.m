  
% &&&&&&&& FUNCTION “scalarfun” &&&&&&&&&&&&&&&
function psi=scalarfun(rx,ry,rz,wk,rad,m,n,del)
    %psi function calculates the vector potential for H-field
    rmx=rx(m);
    rmy=ry(m);
    rmz=rz(m);

    % [rmx,rmy,rmz]=rmhvector(rx,ry,rz,m);
    if del==0.5
        rnx=rx(n);
        rny=ry(n);
        rnz=rz(n);
        [sx1,sy1,sz1,rmag]=sunit(rx,ry,rz,n);
    end
    if del==-0.5
        [rnx,rny,rnz]=rmhvector(rx,ry,rz,n-1);
        [sx1,sy1,sz1,rmag]=sunit(rx,ry,rz,n-1);
    end
    
    delta2=rmag/2.0;
    if m==n
        F=@(s) exp(-j*wk*sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2))./sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2);
        psi=quad(F,0.0,delta2);
    % Self term calculation
    else
    % Non-self term calculation (reduced kernel)
        F=@(s) exp(-j*wk*sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2))./sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2);
        psi=quad(F,0.0,delta2);
    end