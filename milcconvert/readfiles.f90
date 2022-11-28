!Test program to read in gauge field values from milc config files
!Travis H Whyte
!9/15/18

program readfiles

implicit none

real*4, dimension(18,4,24,24,24,64)  :: u
integer :: icri, mu, nx, ny, nz, nt
real*4,  dimension(1) :: beta
real*4, dimension(4) :: tad, alat
integer*4, dimension(4) :: gaction 
integer*4 :: icfgsave, myidin


  open(unit=1,file="/data/whytet/qqcd/milcconvert/l2464f211b600m0102m0509m635a.0002000",action="read",status="old", &
       position="rewind",access="stream")

    read(unit=1) gaction
    read(unit=1) beta
    read(unit=1) alat
    read(unit=1) tad
    read(unit=1) icfgsave
    read(unit=1) myidin

    do icri = 1,18 !changed from 12 to 18 TW 3/16/18
       do mu = 1,4
        do nx = 1,2
         do ny = 1,16
          do nz = 1,24
           do nt = 1,64
           read(unit=1) u(icri,mu,nx,ny,nz,nt)
         enddo ! nt
       enddo ! nz
      enddo !ny
      enddo ! nx
     enddo ! mu 
   enddo ! icri 

    print *, "gaction =", gaction
    print *, "beta =", beta
    print *, "alat =", alat
    print *, "tad =", tad
    print *, "icfgsave =", icfgsave
    print *, "myidin =", myidin
    print *, u(1,1,1,1,1,1)
    close(unit=1,status="keep")

end program readfiles


