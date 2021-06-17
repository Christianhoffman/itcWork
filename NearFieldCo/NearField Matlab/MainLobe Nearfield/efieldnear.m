% &&&&&&&&&&&&& FUNCTION “efield” &&&&&&&&&&&&&&&&&&&&&&&
function efield = efieldnear(rx,ry,rz,xdd,ydd,zdd,solvector,wk,rad,factor,delta,wavelength,ix,iy,iz,nunkns)
    % This function computes the near electric field E_x for ix=1,iy=0,iz=0
    % This function computes the near electric field E_y for ix=0,iy=1,iz=0
    % This function computes the near electric field E_z for ix=0,iy=0,iz=1
    nunkns2=nunkns+2;
    mp1=nunkns2+1;
    if ix == 1
        % **** Compute x-Component of Electric Field (E_x) **********
        rx(mp1)=xdd-0.0005*wavelength;
        ry(mp1)=ydd;
        rz(mp1)=zdd;

        rx(mp1+1)=xdd;
        ry(mp1+1)=ydd;
        rz(mp1+1)=zdd;

        rx(mp1+2)=xdd+0.0005*wavelength;
        ry(mp1+2)=ydd;
        rz(mp1+2)=zdd;
    elseif iy == 1
        rx(mp1)=xdd;
        ry(mp1)=ydd-0.0005*wavelength;
        rz(mp1)=zdd;


        rx(mp1+1)=xdd;
        ry(mp1+1)=ydd;
        rz(mp1+1)=zdd;

        rx(mp1+2)=xdd;
        ry(mp1+2)=ydd+0.0005*wavelength;
        rz(mp1+2)=zdd;
    else
        rx(mp1)=xdd;
        ry(mp1)=ydd;
        rz(mp1)=zdd-0.0005*wavelength;

        rx(mp1+1)=xdd;
        ry(mp1+1)=ydd;
        rz(mp1+1)=zdd;

        rx(mp1+2)=xdd;
        ry(mp1+2)=ydd;
        rz(mp1+2)=zdd+0.0005*wavelength;
    end
    % ****Compute z-Component of Electric Field (E_z) **********

    [rx1,ry1,rz1]= rmhvector(rx,ry,rz,mp1+1);
    [rx2,ry2,rz2]= rmhvector(rx,ry,rz,mp1);

    diffx=rx1-rx2;
    diffy=ry1-ry2;
    diffz=rz1-rz2;

    esum=0.0;

    for n=1:nunkns
        np1=n+1;

        % Contribution due to vector potential
        psi1=vecpot(rx,ry,rz,wk,rad,mp1+1,np1, 0.5);
        psi2=vecpot(rx,ry,rz,wk,rad,mp1+1,np1,-0.5);

        %Contribution due to scalar potential
        psi3=scalarpot(rx,ry,rz,wk,rad,mp1+1,+0.5,np1,np1+1);
        psi4=scalarpot(rx,ry,rz,wk,rad,mp1+1,+0.5,np1-1,np1);
        psi5=scalarpot(rx,ry,rz,wk,rad,mp1+1,-0.5,np1,np1+1);
        psi6=scalarpot(rx,ry,rz,wk,rad,mp1+1,-0.5,np1-1,np1);

        % S unit vectors
        [sx1,sy1,sz1,rmag]=sunit(rx,ry,rz,np1);
        [sx2,sy2,sz2,rmag]=sunit(rx,ry,rz,n);

        % Dot products
        dot1=psi1*(diffx*sx1+diffy*sy1+diffz*sz1);
        dot2=psi2*(diffx*sx2+diffy*sy2+diffz*sz2);
        dotprod=wk^2*(dot1+dot2);

        matrix=factor*(dotprod-psi3/delta+psi4/delta+psi5/delta-psi6/delta);

        esum=esum+matrix*solvector(n);
    end % Loop over n; E_x Component

    efield = -esum/0.001/sqrt(2); % x-comp


