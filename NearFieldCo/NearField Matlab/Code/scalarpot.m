
% &&&&&&&&&& function “scalarpot” &&&&&&&&&&&&&&&&&&&&&&&
function psi=scalarpot(rx,ry,rz,wk,rad,m,del,n,q)
%psi function due to scalar potential
    rnx=rx(n);
    rny=ry(n);
    rnz=rz(n);
    rnx2=rx(q);
    rny2=ry(q);
    rnz2=rz(q);
    [sx1,sy1,sz1,rmag]=sunit(rx,ry,rz,n);
    if del==0.5
        [rmx,rmy,rmz]=rmhvector(rx,ry,rz,m);
    end
    if del==-0.5
        [rmx,rmy,rmz]=rmhvector(rx,ry,rz,m-1);
    end
        delta=rmag;
        delta2=delta/2;
        dist1=sqrt((rmx-rnx)^2+(rmy-rny)^2+(rmz-rnz)^2);
        dist2=sqrt((rmx-rnx2)^2+(rmy-rny2)^2+(rmz-rnz2)^2);
    if (dist1<delta)&&(dist2<delta)
    % Self term calculation
        F=@(s) exp(-j*wk*sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2))./sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2);
        psi=quad(F,0.0,delta);
    else
    % Non-self term calculation - (Reduced kernel Approximation)
        F=@(s) exp(-j*wk*sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2))./sqrt((rmx-rnx-s*sx1).^2+(rmy-rny-s*sy1).^2+(rmz-rnz-s*sz1).^2+rad.^2);
        psi=quad(F,0.0,delta);

    end
