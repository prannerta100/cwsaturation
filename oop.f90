program oop
integer i,a,b,c
open (unit=9,file='baba',status='old',access='sequential')
read(9,*) i
read(9,*) a,b,c

close(unit=9)

open (unit=9,file='baba',status='unknown',access='sequential',form='unformatted')
write(9) 5*i
write(9) a,b,c
close(unit=9)

end program oop

