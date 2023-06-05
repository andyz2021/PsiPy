SUBROUTINE getva(br, br_r, br_t, br_p, bt, bt_r, bt_t, bt_p, bp, bp_r, bp_t, bp_p, rho, rho_r, rho_t, rho_p, nr, nt, np, va)

INTEGER :: i, j, k
INTEGER :: nr, nt, np, br_r, br_t, br_p, bt_r, bt_t, bt_p, bp_r, bp_t, bp_p, rho_r, rho_t, rho_p, nrm1, ntm1
REAL :: br(br_r, br_t, br_p), bt(bt_r, bt_t, bt_p), bp(bp_r, bp_t, bp_p), rho(rho_r, rho_t, rho_p)
REAL, intent(out) :: va(np, nt, nr)
REAL :: avbr, avbt, avbp, bsq, srho


! Calculate v_A
!DO k = 1, np
!    avbr = 0.5 * (br(1, 1, k) + br(1, 1, k))
!    avbt = 0.5 * (bt(1, 1, k) + bt(1, 1, k))
!    avbp = 0.25 * (bp(1, 1, k) + bp(1, 1, k) + bp(1, 1, k) + bp(1, 1, k))
!    bsq = avbr**2 + avbt**2 + avbp**2
!    srho = SIGN(1.0, rho(1, 1, k))
!    va(1, 1, k) = srho * SQRT(bsq / (srho * rho(1, 1, k)))
!END DO
!
!DO k = 1, np
!    DO j = 2, nt
!        avbr = 0.5 * (br(1, j, k) + br(1, j - 1, k))
!        avbt = 0.5 * (bt(1, j, k) + bt(1, j, k))
!        avbp = 0.25 * (bp(1, j, k) + bp(1, j, k) + bp(1, j - 1, k) + bp(1, j - 1, k))
!        bsq = avbr**2 + avbt**2 + avbp**2
!        srho = SIGN(1.0, rho(1, j, k))
!        va(1, j, k) = srho * SQRT(bsq / (srho * rho(1, j, k)))
!    END DO
!END DO
!
!DO k = 1, np
!    DO i = 2, nr
!        avbr = 0.5 * (br(i, 1, k) + br(i, 1, k))
!        avbt = 0.5 * (bt(i, 1, k) + bt(i-1, 1, k))
!        avbp = 0.25 * (bp(i, 1, k) + bp(i-1, 1, k) + bp(i, 1, k) + bp(i, 1, k))
!        bsq = avbr**2 + avbt**2 + avbp**2
!        srho = SIGN(1.0, rho(i, 1, k))
!        va(i, 1, k) = srho * SQRT(bsq / (srho * rho(i, 1, k)))
!    END DO
!END DO
!
!DO k = 1, np
!    DO j = 2, nt
!        DO i = 2, nr
!            avbr = 0.5 * (br(i, j, k) + br(i, j - 1, k))
!            avbt = 0.5 * (bt(i, j, k) + bt(i - 1, j, k))
!            avbp = 0.25 * (bp(i, j, k) + bp(i - 1, j, k) + bp(i, j - 1, k) + bp(i - 1, j - 1, k))
!            bsq = avbr**2 + avbt**2 + avbp**2
!            srho = SIGN(1.0, rho(i, j, k))
!            va(i, j, k) = srho * SQRT(bsq / (srho * rho(i, j, k)))
!        END DO
!    END DO
!END DO

!Original code
nrm1=nr-1
ntm1=nt-1



do i=1,np
    do j=1,nt
        do k=1,nr
            if (k==1 .and. j==1) then
                avbr=br(i,j,k)
                avbt=bt(i,j,k)
                avbp=bp(i,j,k)
            elseif (k==1 .and. j>1) then
                avbr=.5*(br(i,j,k)+br(i,j-1,k))
                avbt=bt(i,j,k)
                avbp = 0.5 * (bp(i, j, k) + bp(i, j - 1, k))
            elseif (j==1 .and. k>1) then
                avbr=br(i,j,k)
                avbt=.5*(bt(i,j,k)+bt(i,j,k-1))
                avbp=0.5*(bp(i,j,k) + bp(i,j,k-1))
            else
                avbr=.5*(br(i,j,k)+br(i,j-1,k))
                avbt=.5*(bt(i,j,k)+bt(i,j,k-1))
                avbp=.25*( bp(i, j, k)+bp(i, j, k-1) + bp(i, j-1, k)+bp(i, j-1, k-1))
            end if
            bsq=(avbr**2+avbt**2+avbp**2)
            srho=sign(1.,rho(i,j,k))
            va(i,j,k)=srho*sqrt(bsq/(srho*rho(i,j,k)))
        enddo
    enddo
enddo

END SUBROUTINE getva