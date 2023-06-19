! ============================================================================
! Name        : readfort10ij.f90
! Author      : Biaobin Li
! Version     : V-1.0
! Copyright   : Your copyright notice
! Description : combine fort.10ij into one file:
!               impt_phase.txt, impt_id.txt
! ============================================================================

program readfort10ij
    implicit none
    integer :: j,fortstart,fortend,np
    integer :: fort,stat,cnt
    character(len=100) :: filename
    real*8 :: x,bgx,y,bgy,z,bgz

    !==============================
    print*,"Enter the start and end file number:"
    read*,fortstart,fortend
    !fortstart = 1000
    !fortend   = 1089
    !==============================


    np = fortend-fortstart+1

    open(10,file="impt_phase.txt")
    open(11,file="impt_id.txt")

    cnt = 1
    write(11,*)cnt
    do j=1,np
        fort = fortstart+j-1
        write(filename,*)fort
        filename = "fort."//adjustl(filename)

        print*,filename

        open(2,file=filename)
        do while (.true.)
            read(2,*,iostat=stat)x,bgx,y,bgy,z,bgz

            !exit the file when at the end
            if(stat/=0) exit

            !conver to cm
            x = x*100.0
            y = y*100.0
            z = z*100.0
            write(10,101)x,bgx,y,bgy,z,bgz

            cnt = cnt+1
        enddo
        print*,"cnt=",cnt
        write(11,*)cnt
    enddo

101 format(6(1x,e20.12))

close(2)
close(10)
close(11)

print*,"data combine finished."

end program
