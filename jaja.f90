program jaja
    integer :: ndimo, i
    character(100) ixname
    integer, allocatable, dimension(:) :: mpi1,mqi1,ml1,mjk1,mk1,mjm1,mm1
    if(command_argument_count().ne.1) stop
    call get_command_argument(1,ixname)
    !write(*,*)'opening ixname ', ixname
    open (unit=9,file=ixname,status='old', access='sequential',form='unformatted')
    !write(*,*)'opened it'
    read (9) ndimo
    allocate(mpi1(ndimo),mqi1(ndimo),ml1(ndimo),mjk1(ndimo),mk1(ndimo),mjm1(ndimo),mm1(ndimo))
    read (9) (mpi1(i),mqi1(i),ml1(i),mjk1(i),mk1(i),mjm1(i),mm1(i),i=1,ndimo)
    close (unit=9)
    write(*,*) ndimo
    do i=1,ndimo
        write(*,*) mpi1(i),mqi1(i),ml1(i),mjk1(i),mk1(i),mjm1(i),mm1(i)
    enddo
end

