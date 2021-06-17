% z-comp of Electric Field (Near Field Calculation for a Wire
%Antenna) for near field calculation August, 2007
clear all
%******************** INPUT *******************************                 % Variables with asterik (*) on the comment on the right side are INPUT
nunkns=31;                                                                  % ******** Number of Unknowns (for current) on the Antenna
nunkns2=nunkns+2;
freq=2.4*10^9;                                                              % ******************* Frequency (Hz)
vel=3*10^8;                                                                 % velocity of light in free space
omega=2*pi*freq;                                                            % angular frequency (rad/sec)
wk=omega/vel;                                                               % propagation constant (2*pi/lambda)
wavelength=vel/freq;                                                        % wavelength (m)
length=0.2242*wavelength;                                                   % ****************** Half dipole length (in lambda)
dlength=2*length;
rad=0.005*wavelength;                                                       % ******************** Radius of Dipole (inlambda)
waveimp=377;                                                                % ********** Wave Impedance in free space(ohms)
eps=1/(36*pi*10^9);                                                         % ********** Free space permittivity (F/m)
zg=0;                                                                       % ********* Feed point of the antenna (zg=0 for Center-fed)
delta=2*length/(nunkns+1);                                                  % Antenna segment length (2*L/(N+1))
                                                                            % Grid points (x,y,z) where near field needs to be computed
dx=0;
dy=0;
dz=0.001;
ndx=1;
ndy=1;
ndz=300;
xstart=0.0*wavelength;
ystart_array(1)=0.01;
ystart_array(2)=0.03;
ystart_array(3)=0.05;
zstart=0.001*wavelength;
                                                                            % ****** For Antenna case, the Dipole is Center-fed with 1 Volt ******
                                                                            % Forcing Function "vmvector" is Unity at feedpoint
for p=1:nunkns
    matchpnt=((nunkns-1)/2)+1;
    if p==matchpnt
        vmvector(p,1)=10; %12.25;                                                          % ******* Antenna feed voltage is 1 volt
    else
    vmvector(p,1)=0;
    end
end

% ************************** End of Input Data********************
% The array rx,ry,rz gives the x,y,z components of wire segment end points
for nn=1:nunkns2
    rx(nn)=0.0;
    ry(nn)=0.0;
    rz(nn)=-length+delta*(nn-1);
end
% Calculation of the Impedance Matrix "zmatrix"
factor=-1/(j*4*pi*omega*eps);
for m=1:nunkns
    mp1=m+1;
    [rx1,ry1,rz1]= rmhvector(rx,ry,rz,mp1);
    [rx2,ry2,rz2]= rmhvector(rx,ry,rz,m);
    diffx=rx1-rx2;
    diffy=ry1-ry2;
    diffz=rz1-rz2;
    for n=1:nunkns
        np1=n+1;

        % Contribution due to vector potential

        psi1=vecpot(rx,ry,rz,wk,rad,mp1,np1, 0.5);
        psi2=vecpot(rx,ry,rz,wk,rad,mp1,np1,-0.5);

        %Contribution due to scalar potential

        psi3=scalarpot(rx,ry,rz,wk,rad,mp1,+0.5,np1,np1+1);
        psi4=scalarpot(rx,ry,rz,wk,rad,mp1,+0.5,np1-1,np1);
        psi5=scalarpot(rx,ry,rz,wk,rad,mp1,-0.5,np1,np1+1);
        psi6=scalarpot(rx,ry,rz,wk,rad,mp1,-0.5,np1-1,np1);

        % S unit vectors
        [sx1,sy1,sz1,rmag]=sunit(rx,ry,rz,np1);
        [sx2,sy2,sz2,rmag]=sunit(rx,ry,rz,n);

        % Dot products
        dot1=psi1*(diffx*sx1+diffy*sy1+diffz*sz1);
        dot2=psi2*(diffx*sx2+diffy*sy2+diffz*sz2);
        dotprod=wk^2*(dot1+dot2);


        zmatrix(m,n)=factor*(dotprod-psi3/delta+psi4/delta+psi5/delta-psi6/delta);

    end
end

% SOLUTION BY MATRIX INVERSION

solvector=inv(zmatrix)*vmvector;                                            % solvector is the solution vector

%with current on the antenna
% Plot the Current Distribution on the Wire Antenna

rsolvec=real(solvector);
isolvec=imag(solvector);
for ip=1:nunkns
    rrealpart(ip)=rsolvec(ip);
    iimagpart(ip)=isolvec(ip);
    zpart(ip)=rz(ip);
end

for ipp=1:nunkns2
    if ipp==1
        realpart(ipp)=0;
        imagpart(ipp)=0;
    elseif ipp==nunkns2
        realpart(ipp)=0;
        imagpart(ipp)=0;
    else
        realpart(ipp)=rrealpart(ipp-1);
        imagpart(ipp)=iimagpart(ipp-1);
    end
end


% Near Field Computations
mp1=nunkns2+1;
for kk=1:3
    nfield_points=0;
    ystart=ystart_array(kk);
    for i=1:ndx
        xd(i)=xstart+(i-1)*dx;
        xdi=xd(i);
        for j=1:ndy
            yd(j)=ystart+(j-1)*dy;
            ydj=yd(j);
            for k=1:ndz
                nfield_points=nfield_points+1;
                zd(k)=zstart+(k-1)*dz;
                zdk=zd(k);
                % efieldx,efieldy and efieldz are the x,y,z components of the
                %Electric Field

                % **** NOTE ***** In this example only the z-comp of E-field is
                %being computed
                % therefore, efieldx and efieldy are commented out

                %efieldx=efieldnear(rx,ry,rz,xdi,ydj,zdk,solvector,wk,rad,factor,delta,
                %wavelength,1,0,0,nunkns);

                %efieldy=efieldnear(rx,ry,rz,xdi,ydj,zdk,solvector,wk,rad,factor,delta,
                %wavelength,0,1,0,nunkns);

                efieldz=efieldnear(rx,ry,rz,xdi,ydj,zdk,solvector,wk,rad,factor,delta,wavelength,0,0,1,nunkns);

                % nearfield(kk,nfield_points)=abs(efieldx); % if x-comp of near
                %field is needed
                % nearfield(kk,nfield_points)=abs(efieldy); % if y-comp of near
                %field is needed
                nearfield(kk,nfield_points)=abs(efieldz); % if z-comp of near
                %field is needed

            end % loop over k (z values)
        end % loop for j (y values)
    end % loop over i (x values)
end % loop over kk



% If following three lines are commented out The Antenna Current is
%not plotted
 %plot(rz,realpart,'k-',rz,imagpart,'o')
 %grid on
 %title('Current Distribution over a lamda/2 dipole @ 2.4GHz')
 %xlabel('Length (m)')
 %ylabel('Current Distribution (A)')



 Plot of Electric Field vs z/lambda for three values of rho/lambda
% % % 
plot(zd(:),nearfield(1,:),'k-')
hold on
grid on
%plot(zd(:),nearfield(2,:),'b-')

%plot(zd(:),nearfield(3,:),'r-')
title('Z-component of Electric Field for lamda/2 dipole @.......Hz')
xlabel('Z/lambda')
ylabel('Electric Field in (V/m)')
%gtext('y (lambda) =0.01')
%gtext('0.03')
%gtext('0.05')
hold off