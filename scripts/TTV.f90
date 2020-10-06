implicit real*8(a-h,o-z)

character*50 file001,file002
character*50 j1,j2,j3,j4,j5,j6,j7,j8,j9,j10

file001="case5_history.txt"
file002="case5_transit_times.txt"

open(1,file=file001,status="unknown")
open(2,file=file002,status="unknown")

iflag_transit=0

read(1,*)j1,j2,j3,j4,j5,j6,j7,j8,j9,j10

nn=0
i=1
12	  read(1,*,end=800)z,iii,z,z,z,z,z,z,z
        nn=nn+1
        i=i+1
        goto 12
 800    continue
write(*,*)nn
close(1)

open(1,file=file001,status="unknown")

read(1,*)j1,j2,j3,j4,j5,j6,j7,j8,j9,j10

do i=1,nn
    read(1,*)ti,iii,zero,xi,zero,zero,zero,zero,zero
    if(i.gt.(1)) then
        if(xi.gt.(0).and.xim1.lt.(0)) then
            zl1=-xi/(xim1-xi)
            zl2=-xim1/(xi-xim1)
            time_transit=tim1*zl1+ti*zl2
            write(2,*)iflag_transit,time_transit
            iflag_transit=iflag_transit+1
        endif
    endif
    tim1=ti
    xim1=xi
enddo

end
