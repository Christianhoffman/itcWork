% &&&&&&&&&&&&& FUNCTION “sunit” &&&&&&&&&&&&&&&&&&
function [snx,sny,snz,rmag] = sunit(rx,ry,rz,p)
    %Gives x,y,z components of the unit vector,s_(p+1/2), on segment p to
    %p+1
    rx1=rx(p+1);
    ry1=ry(p+1);
    rz1=rz(p+1);
    rx2=rx(p);
    ry2=ry(p);
    rz2=rz(p);
    rmag=sqrt((rx1-rx2)^2+(ry1-ry2)^2+(rz1-rz2)^2);
    snx=(rx1-rx2)/rmag;
    sny=(ry1-ry2)/rmag;
    snz=(rz1-rz2)/rmag;
    