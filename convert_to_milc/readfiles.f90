program readfiles

implicit none

real*4, dimension(18,192,4,2,16)  :: u
integer :: icri, isite, mu, ieo, ibl
real*4,  dimension(1) :: beta
real*4, dimension(4) :: tad, alat
integer*4, dimension(4) :: gaction 
integer*4 :: icfgsave, myidin


  open(unit=1,file="l2464f211b600m0102m0509m635a.0001004",action="read",status="old", &
       position="rewind",form="binary",access="stream")

    read(unit=1) gaction
    read(unit=1) beta
    read(unit=1) alat
    read(unit=1) tad
    read(unit=1) icfgsave
    read(unit=1) myidin

    do icri = 1,18 !changed from 12 to 18 TW 3/16/18
      do isite = 1,192
       do mu = 1,4
        !do ieo = 1,2
         do ibl = 1,16
        read(unit=1) (u(icri,isite,mu,1,ibl), u(icri,isite,mu,2,ibl))
           print *, u(icri,isite,mu,1,ibl), u(icri,isite,mu,2,ibl)
       enddo ! ibl
       !enddo ! ieo
      enddo !mu
      enddo ! isite
     enddo ! icri  

    print *, "gaction =", gaction
    print *, "beta =", beta
    print *, "alat =", alat
    print *, "tad =", tad
    print *, "icfgsave =", icfgsave
    print *, "myidin =", myidin
    close(unit=1,status="keep")

end program readfiles


