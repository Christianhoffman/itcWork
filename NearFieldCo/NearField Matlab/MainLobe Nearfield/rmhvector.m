% &&&&&&&&&& FUNCTION “rmhvector” &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function[rmhx,rmhy,rmhz] = rmhvector(rx,ry,rz,p)
    %Gives the x,y,z components of vector locating the mid-point of a
    %segment
    rmhx=(rx(p+1)+rx(p))/2;
    rmhy=(ry(p+1)+ry(p))/2;
    rmhz=(rz(p+1)+rz(p))/2;