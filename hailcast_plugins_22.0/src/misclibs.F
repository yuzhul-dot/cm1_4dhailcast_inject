  MODULE misclibs

  implicit none

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setrho(rho,rr,rf,c1,c2,                                 &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                 &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_s)

      use input, only : ib,ie,jb,je,kb,ke,cmp,imp,jmp,kmp,kmt,rmp,ni,nj,nk,ngt, &
                        cgs1,cgs2,cgs3,cgt1,cgt2,cgt3
      use bc_module, only : bcs
#ifdef MPI
      use bc_module, only : bcs2
      use comm_module, only : comm_1s_start,comm_1s_end,getcorner
#endif
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(inout), dimension(rmp) :: reqs_s

      integer :: i,j,k

      call bcs(rho)
#ifdef MPI
      call comm_1s_start(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
      call comm_1s_end(  rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
      call getcorner(    rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcs2(rho)
#endif

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=2,nk
      do j=0,nj+1
      do i=0,ni+1
        rr(i,j,k) = 1.0/rho(i,j,k)
        rf(i,j,k) = (c1(i,j,k)*rho(i,j,k-1)+c2(i,j,k)*rho(i,j,k))
      enddo
      enddo
      enddo

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        rr(i,j,1) = 1.0/rho(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        rf(i,j,1) = cgs1*rho(i,j,1)+cgs2*rho(i,j,2)+cgs3*rho(i,j,3)
        rf(i,j,nk+1) = cgt1*rho(i,j,nk)+cgt2*rho(i,j,nk-1)+cgt3*rho(i,j,nk-2)
      enddo
      enddo

      end subroutine setrho


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


      subroutine setuvw(u3d,v3d,w3d,n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,setw)

      use input, only : ib,ie,jb,je,kb,ke,cmp,kmt
      use bc_module, only : bcu2,bcv2,bcw2
      use comm_module, only : getcorneru3,getcornerv3,getcornerw3
      implicit none

      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      integer, intent(in) :: setw

      call getcorneru3(u3d,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                           s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call getcornerv3(v3d,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                           s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      if( setw.eq.1 ) call getcornerw3(w3d,n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2)

      call bcu2(u3d)
      call bcv2(v3d)
      if( setw.eq.1 ) call bcw2(w3d)

      end subroutine setuvw


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tavg_plus1(rdenom,tavg,timavg,var,istag,jstag,kstag,k1,k2,dbldt)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,qp,ibta,ieta,jbta,jeta,ntim

      implicit none

      ! add one more data point to time-averages

      double precision, intent(in) :: rdenom
      integer, intent(in) :: istag,jstag,kstag,k1,k2
      real,             intent(inout), dimension(ibta:ieta,jbta:jeta,k1:k2,ntim) :: tavg
      real,             intent(inout), dimension(ibta:ieta,jbta:jeta,k1:k2) :: timavg
      real, intent(in), dimension(ib:ie+istag,jb:je+jstag,k1:k2+kstag) :: var
      double precision, intent(in) :: dbldt

      integer :: i,j,k,n
      real(kind=qp) :: qtem

      IF( ntim.eq.5 )THEN

        do k=k1,k2
        do j=-1,nj+2
        do i=-1,ni+2
          tavg(i,j,k,ntim) = tavg(i,j,k,ntim) + var(i,j,k)*dbldt
          timavg(i,j,k) = rdenom*( tavg(i,j,k,1) &
                                  +tavg(i,j,k,2) &
                                  +tavg(i,j,k,3) &
                                  +tavg(i,j,k,4) &
                                  +tavg(i,j,k,5) )
        enddo
        enddo
        enddo

      ELSE

        do k=k1,k2
        do j=-1,nj+2
        do i=-1,ni+2
          tavg(i,j,k,ntim) = tavg(i,j,k,ntim) + var(i,j,k)*dbldt
          qtem = 0.0
          do n=1,ntim
            qtem = qtem + tavg(i,j,k,n)
          enddo
          timavg(i,j,k) = rdenom*qtem
        enddo
        enddo
        enddo

      ENDIF

      end subroutine tavg_plus1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tavg_shift(tavg,k1,k2)

      use input, only : ni,nj,ibta,ieta,jbta,jeta,ntim

      implicit none

      ! add one more data point to time-averages

      integer, intent(in) :: k1,k2
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,k1:k2,ntim) :: tavg

      integer :: i,j,k,n

      DO n = 1 , ntim-1

        do k=k1,k2
        do j=-1,nj+2
        do i=-1,ni+2
          tavg(i,j,k,n) = tavg(i,j,k,n+1)
        enddo
        enddo
        enddo

      ENDDO

      do k=k1,k2
      do j=-1,nj+2
      do i=-1,ni+2
        tavg(i,j,k,ntim) = 0.0
      enddo
      enddo
      enddo

      end subroutine tavg_shift


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getdiv(arh1,arh2,uh,vh,mh,u,v,w,dum1,dum2,dum3,div,  &
                        rds,rdsf,sigma,sigmaf,gz,rgzu,rgzv,dzdx,dzdy)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,itb,ite,jtb,jte,  &
                        terrain_flag,axisymm,                        &
                        rdx,rdy,rdz,smeps,zt,rzt,                    &
                        timestats,time_divx,mytime

      implicit none

      real, intent(in), dimension(ib:ie) :: arh1,arh2,uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,div
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgzu,rgzv,dzdx,dzdy

      integer :: i,j,k
      real :: r1,r2

      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          ! Cartesian without terrain:
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            div(i,j,k)=( (u(i+1,j,k)-u(i,j,k))*rdx*uh(i)        &
                        +(v(i,j+1,k)-v(i,j,k))*rdy*vh(j) )      &
                        +(w(i,j,k+1)-w(i,j,k))*rdz*mh(1,1,k)
            if(abs(div(i,j,k)).lt.smeps) div(i,j,k)=0.0
          enddo
          enddo
          enddo
        ELSE
          ! axisymmetric:
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            div(i,j,k)=(arh2(i)*u(i+1,j,k)-arh1(i)*u(i,j,k))*rdx*uh(i)   &
                      +(w(i,j,k+1)-w(i,j,k))*rdz*mh(1,1,k)
            if(abs(div(i,j,k)).lt.smeps) div(i,j,k)=0.0
          enddo
          enddo
          enddo
        ENDIF
      ELSE
          ! Cartesian with terrain:
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          DO k=1,nk
            do j=1,nj
            do i=1,ni+1
              dum1(i,j,k)=u(i,j,k)*rgzu(i,j)
            enddo
            enddo
            do j=1,nj+1
            do i=1,ni
              dum2(i,j,k)=v(i,j,k)*rgzv(i,j)
            enddo
            enddo
          ENDDO
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k,r1,r2)
          DO k=1,nk
            IF(k.eq.1)THEN
              do j=1,nj
              do i=1,ni
                dum3(i,j,1)=0.0
                dum3(i,j,nk+1)=0.0
              enddo
              enddo
            ELSE
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              r1 = 0.5*r1
              r2 = 0.5*r2
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=w(i,j,k)                                             &
                           +( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                               )*(sigmaf(k)-zt)*gz(i,j)*rzt
              enddo
              enddo
            ENDIF
          ENDDO
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            div(i,j,k)=( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)        &
                        +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j) )      &
                        +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k)
            if(abs(div(i,j,k)).lt.smeps) div(i,j,k)=0.0
          enddo
          enddo
          enddo
      ENDIF
      if(timestats.ge.1) time_divx=time_divx+mytime()

      end subroutine getdiv


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getdivx(arh1,arh2,uh,vh,mh,rho0,rf0,rru,rrv,rrw,divx,  &
                         rds,rdsf,sigma,sigmaf,gz,rgzu,rgzv,dzdx,dzdy)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,itb,ite,jtb,jte,  &
                        terrain_flag,axisymm,                        &
                        rdx,rdy,rdz,smeps,zt,rzt,                    &
                        timestats,time_divx,time_advs,mytime

      implicit none

      real, intent(in), dimension(ib:ie) :: arh1,arh2,uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rf0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgzu,rgzv,dzdx,dzdy

      integer :: i,j,k
      real :: r1,r2,tem

    IF(.not.terrain_flag)THEN
      ! without terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
      DO k=1,nk
        tem = rho0(1,1,k)
        do j=1,nj
        do i=1,ni+1
          rru(i,j,k)=rru(i,j,k)*tem
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          rrv(i,j,k)=rrv(i,j,k)*tem
        enddo
        enddo
        IF(k.eq.1)THEN
          do j=1,nj
          do i=1,ni
            rrw(i,j,   1) = 0.0
            rrw(i,j,nk+1) = 0.0
          enddo
          enddo
        ELSE
          tem = rf0(1,1,k)
          do j=1,nj
          do i=1,ni
            rrw(i,j,k)=rrw(i,j,k)*tem
          enddo
          enddo
        ENDIF
      ENDDO

    ELSE
      ! with terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      DO k=1,nk
        do j=1,nj
        do i=1,ni+1
          rru(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*rru(i,j,k)*rgzu(i,j)
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          rrv(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*rrv(i,j,k)*rgzv(i,j)
        enddo
        enddo
      ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,r2)
      DO k=1,nk
        IF(k.eq.1)THEN
          do j=1,nj
          do i=1,ni
            rrw(i,j,   1) = 0.0
            rrw(i,j,nk+1) = 0.0
          enddo
          enddo
        ELSE
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          r1 = 0.5*r1
          r2 = 0.5*r2
          do j=1,nj
          do i=1,ni
            rrw(i,j,k)=rf0(i,j,k)*rrw(i,j,k)                              &
                      +( ( r2*(rru(i,j,k  )+rru(i+1,j,k  ))               &
                          +r1*(rru(i,j,k-1)+rru(i+1,j,k-1)) )*dzdx(i,j)   &
                        +( r2*(rrv(i,j,k  )+rrv(i,j+1,k  ))               &
                          +r1*(rrv(i,j,k-1)+rrv(i,j+1,k-1)) )*dzdy(i,j)   &
                       )*(sigmaf(k)-zt)*gz(i,j)*rzt
          enddo
          enddo
        ENDIF
      ENDDO

    ENDIF
    if(timestats.ge.1) time_advs=time_advs+mytime()

      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          ! Cartesian without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                         +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )      &
                         +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ELSE
          ! axisymmetric:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)   &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
        ENDIF
      ELSE
          ! Cartesian with terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k)=( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                         +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )      &
                         +(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo
          enddo
      ENDIF
      if(timestats.ge.1) time_divx=time_divx+mytime()

      end subroutine getdivx


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine randpert(xfref,yfref,xh,yh,zh,zf,tha)

      ! Generate random perturbations, add to tha (by default):

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,nx,ny,ngxy

      implicit none

      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha

      integer :: i,j,k
      real :: rand,amplitude


      ! reinitialize pseudorandom number generator:
      call reinit_random_seed


      ! this is the amplitude of the theta perturbations
      ! (plus or minus this value in K)
      amplitude = 0.25


      ! random numbers added here
      ! (can be modified to only place perturbations in certain
      !  locations, but this default code simply puts them
      !  everywhere)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        call random_number(rand)
        tha(i,j,k) = tha(i,j,k)+amplitude*(2.0*rand-1.0)
      enddo
      enddo
      enddo

      end subroutine randpert


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine reinit_random_seed

      use input, only : myid

      implicit none

      integer :: k,n
      integer,dimension(8) :: values
      integer, dimension(:), allocatable :: sand
      real :: rand

      call date_and_time(values=values)
      call random_seed(size=k)
      k = max(3,k)
      allocate( sand(k) )
      sand(1) = myid
      sand(2) = values(7)
      sand(3) = values(8)
      call random_seed(put=sand(1:k))
      call random_number(rand)
      print *,'  myid,rand = ',myid,rand
      deallocate( sand )

      end subroutine reinit_random_seed


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)

      implicit none

      integer, intent(in) :: myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe
      real, intent(in) :: zdeep,lamx,lamy,xcent,ycent,aconv
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in),    dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d

      integer :: i,j,k
      real :: term1,term2,term3,term4,umo

!!!      if(myid.eq.0) print *,'    convinitu '
!$omp parallel do default(shared)   &
!$omp private(i,j,k,term1,term2,term3,term4,umo)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        term4 = (zdeep-0.5*(zh(i-1,j,k)+zh(i,j,k)))/zdeep
        if (term4 .gt. 0.0) then
          term1 = -(2.0*Aconv*(xf(i)-xcent))/(lamx**2)
          term2 = -((xf(i)-xcent)/lamx)**2
          term3 = -((yh(j)-ycent)/lamy)**2
          umo = term1*(exp(term2)*exp(term3))*term4
          if( abs(umo).gt.0.01 ) u3d(i,j,k) = u0(i,j,k)+umo
        endif
      enddo
      enddo
      enddo

      end subroutine convinitu


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
      implicit none

      integer, intent(in) :: myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn
      real, intent(in) :: zdeep,lamx,lamy,xcent,ycent,aconv
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in),    dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d

      integer :: i,j,k
      real :: term1,term2,term3,term4,vmo

!!!      if(myid.eq.0) print *,'    convinitv '
!$omp parallel do default(shared)   &
!$omp private(i,j,k,term1,term2,term3,term4,vmo)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        term4 = (zdeep-0.5*(zh(i,j-1,k)+zh(i,j,k)))/zdeep
        if (term4 .gt. 0.0) then
          term1 = -(2.0*Aconv*(yf(j)-ycent))/(lamy**2)
          term2 = -((xh(i)-xcent)/lamx)**2
          term3 = -((yf(j)-ycent)/lamy)**2
          vmo = term1*(exp(term2)*exp(term3))*term4
          if( abs(vmo).gt.0.01 ) v3d(i,j,k) = v0(i,j,k)+vmo
        endif
      enddo
      enddo
      enddo

      end subroutine convinitv


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine get_wnudge(mtime,dtin,xh,yh,zf,wa,wndgten)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,  &
                        t1_wnudge,t2_wnudge,xc_wnudge,yc_wnudge,zc_wnudge,alpha_wnudge,  &
                        wmax_wnudge,rxrwnudge,ryrwnudge,rzrwnudge
      use constants , only : pi

      implicit none

      double precision, intent(in) :: mtime
      real, intent(in) :: dtin
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: wndgten

      integer :: i,j,k
      real :: beta,wmag,gamm,tem

      !  updraft nudging scheme (Naylor and Gilmore, 2012, MWR, pgs 3699-3705)

      gamm = 1.0

      if(mtime.ge.t1_wnudge)THEN
        gamm = 1.0+(0.0-1.0)*(mtime-t1_wnudge)/(t2_wnudge-t1_wnudge)
      endif

!!!      if(myid.eq.0) print *,'    get_wnudge: mtime,gamm = ',mtime,gamm

      tem = alpha_wnudge * gamm

!$omp parallel do default(shared)   &
!$omp private(i,j,k,beta,wmag)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        beta = sqrt( ((xh(i)-xc_wnudge)*rxrwnudge)**2       &
                    +((yh(j)-yc_wnudge)*ryrwnudge)**2       &
                    +((zf(i,j,k)-zc_wnudge)*rzrwnudge)**2)
        if(beta.lt.1.0)then
          wmag = wmax_wnudge*( cos(0.5*pi*beta)**2 )
          wndgten(i,j,k) = tem*max(wmag-wa(i,j,k),0.0)
        else
          wndgten(i,j,k) = 0.0
        endif
      enddo
      enddo
      enddo

      end subroutine get_wnudge


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getnewdt(ndt,dt,dtlast,adt,acfl,dbldt,                                 &
                          mtime,stattim,taptim,rsttim,prcltim,hailtrajtim,     &
                                                  diagtim,azimavgtim,       &
                          dorestart,dowriteout,dostat,doprclout,dohailout,dotdwrite,doazimwrite,  &
                          hifrqtim,dohifrqwrite,doinit)

      use input, only : nx,ny,dtl,ksmax,cflmax,cm1setup,imoist,ptype,tapfrq,statfrq,rstfrq, &
                        myid,dodomaindiag,diagfrq,doazimavg,azimavgfrq,dohifrq,hifrqfrq,    &
                        smeps,cfl_limit,ramp_down,timestats,mytime,time_misc,time_microphy
      use goddard_module, only : consat2
      use lfoice_module, only : lfoice_init
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(inout) :: ndt
      real, intent(inout) :: dt,dtlast
      double precision, intent(inout) :: adt,acfl,dbldt
      double precision, intent(in) :: mtime,stattim,taptim,rsttim,prcltim,diagtim,azimavgtim,hailtrajtim
      logical, intent(in) :: dorestart,dowriteout,dostat,doprclout,dotdwrite,doazimwrite,dohailout
      double precision, intent(in) :: hifrqtim
      logical, intent(in) :: dohifrqwrite
      logical, intent(in) :: doinit

      integer :: ierr
      real :: tem,ks_limit
      double precision :: tout

             ! cm1r22.0: this is now a namelist parameter 
!!!      real, parameter  ::  cfl_limit   =  1.00    ! maximum CFL allowed  (actually a "target" value)

      real, parameter  ::  max_change  =  0.10    ! maximum (percentage) change in timestep

      if( cm1setup.ge.1 )then
        ! assume vertical component is handled implicitly
        if(nx.gt.3.and.ny.gt.3)then
          ! 3d:
          ks_limit = 0.16/2.0
        else
          ! 2d (including axisymm):
          ks_limit = 0.18
        endif
      endif

      dtlast = dt

      myid0:  &
      IF( myid.eq.0 )THEN
        ! only processor 0 does this:

        cflmax = max(cflmax,1.0e-10)

        IF( cflmax.gt.cfl_limit )THEN
          ! decrease timestep:
          dbldt = dbldt*(cfl_limit/cflmax)
        ELSE
          ! increase timestep:
          dbldt = dbldt*min( 1.0+max_change , cfl_limit/cflmax )
        ENDIF

        ! 180129:
        IF( cm1setup.ge.1 )THEN
        if( (ksmax*dbldt/dt).gt.ks_limit )then
          ! decrease timestep:
          dbldt = min( dbldt , dble(dt*ks_limit/max(1.0e-10,ksmax)) )
        endif
        ENDIF

        ! don't allow dt to exceed twice initial timestep
        dbldt = min( dbldt , dble(2.0*dtl) )


        IF( (ramp_down.eq.1) .and. (.not. doinit) )THEN

            IF( taptim.gt.0.0 .and. tapfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching output time
              if( dowriteout )then
                tout = ( ( taptim - mtime ) + tapfrq )
              else
                tout = ( taptim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming writeout time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001d0*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming writeout time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming writeout time '
                dbldt = tout
              endif
            ENDIF

            IF( rsttim.gt.0.0 .and. rstfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching restart time
              if( dorestart )then
                tout = ( ( rsttim - mtime ) + rstfrq )
              else
                tout = ( rsttim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming restart time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001d0*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming restart time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming restart time '
                dbldt = tout
              endif
            ENDIF

            IF( stattim.gt.0.0 .and. statfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching stat time
              if( dostat )then
                tout = ( ( stattim - mtime ) + statfrq )
              else
                tout = ( stattim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming stat time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001d0*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming stat time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming stat time '
                dbldt = tout
              endif
            ENDIF

            IF( dodomaindiag .and. diagtim.gt.0.0 .and. diagfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching domaindiag time
              if( dotdwrite )then
                tout = ( ( diagtim - mtime ) + diagfrq )
              else
                tout = ( diagtim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming diag time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001d0*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming diag time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming diag time '
                dbldt = tout
              endif
            ENDIF

            IF( doazimavg .and. azimavgtim.gt.0.0 .and. azimavgfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching azimavg time
              if( doazimwrite )then
                tout = ( ( azimavgtim - mtime ) + azimavgfrq )
              else
                tout = ( azimavgtim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming azimavg time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming azimavg time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming azimavg time '
                dbldt = tout
              endif
            ENDIF

            IF( dohifrq .and. hifrqtim.gt.0.0 .and. hifrqfrq.gt.0.0 )THEN
              ! ramp-down timestep when approaching hifrq time
              if( dohifrqwrite )then
                tout = ( ( hifrqtim - mtime ) + hifrqfrq )
              else
                tout = ( hifrqtim - mtime )
              endif
              if( tout.gt.(2.00001d0*dbldt) .and. tout.le.(3.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming hifrq time '
                dbldt = tout/3.0d0
              elseif( tout.gt.(1.00001d0*dbldt) .and. tout.le.(2.0d0*dbldt)  )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming hifrq time '
                dbldt = tout/2.0d0
              elseif( tout.le.dbldt )then
                if(myid.eq.0) print *,'  Note: setting dt based on upcoming hifrq time '
                dbldt = tout
              endif
            ENDIF

        ENDIF

        ! end of processor 0 stuff
      ENDIF  myid0


      ! all processors:
#ifdef MPI
      call MPI_BCAST(dbldt ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

      dt = dbldt
      call dtsmall(dt,dbldt)

      if( dt.lt.smeps )then
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        if(myid.eq.0)then
          print *
          print *,'  getnewdt: dt,dbldt = ',dt,dbldt
          print *
          print *,'  stopping cm1  '
          print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

      IF( ( dt.ne.dtlast ) .or. doinit )THEN
        IF( (imoist.eq.1).and.(ptype.eq.2) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call consat2(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF
        IF( (imoist.eq.1).and.(ptype.eq.4) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call lfoice_init(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF
      ENDIF

      tem = dt/dtlast
      cflmax = cflmax*tem
      ksmax = ksmax*tem

      ndt = ndt + 1
      adt = adt + dbldt
      acfl = acfl + cflmax

      if(timestats.ge.1) time_misc=time_misc+mytime()

      end subroutine getnewdt


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine dtsmall(dt,dbldt)

      use input, only : myid,adapt_dt,psolver,nsound,nx,ny,min_dx,min_dy,min_dz,csound
!      use constants, only : 
#ifdef MPI
      use mpi
#endif
      implicit none

      ! cm1r18:  moved this section of code from solve.F to misclibs.F

      real, intent(inout) :: dt
      double precision, intent(inout) :: dbldt

      integer :: ierr
      real :: dtsm

      ! GHB:  this value is arbitrary ... could be probably be changed
!!!      integer, parameter :: max_nsound = 96
      integer, parameter :: max_nsound = 192

      IF( psolver.eq.1 .or. psolver.eq.4 .or. psolver.eq.5 )THEN

        dtsm = dt
        nsound = 1

      ELSE

        ! Algorithm to determine number of small steps:
        IF( psolver.eq.2 )THEN
          ! check dx,dy,dz:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.60*min( min_dx , min_dz )/360.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.60*min( min_dy , min_dz )/360.0
          ELSE
            ! 3D sims:
            dtsm = 0.50*min( min_dx , min_dy , min_dz )/360.0
          ENDIF
        ELSEIF( psolver.eq.3 )THEN
          ! check dx,dy:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.60*min_dx/360.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.60*min_dy/360.0
          ELSE
            ! 3D sims:
            dtsm = 0.60*min( min_dx , min_dy )/360.0
          ENDIF
        ELSEIF( psolver.eq.6 .or. psolver.eq.7 )THEN
          ! check dx,dy,dz:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.60*min( min_dx , min_dz )/csound
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.60*min( min_dy , min_dz )/csound
          ELSE
            ! 3D sims:
            dtsm = 0.50*min( min_dx , min_dy , min_dz )/csound
          ENDIF
        ENDIF

        nsound = max( nint( dbldt/dtsm ) , 4 )
        if( mod(nsound,2).ne.0 ) nsound = nsound + 1
        if( dbldt/float(nsound).gt.dtsm ) nsound = nsound + 2

        if( nsound.gt.max_nsound )then
          ! GHB:  this is arbitrary ... could be changed
          if( adapt_dt.eq.1 )then
            nsound = max_nsound
            dbldt = nsound*dtsm
            dt = dbldt
          else
            if(myid.eq.0)then
            print *,'  -------------------------------- '
            print *
            print *,'  Limit for number of small steps exceeded: '
            print *
            print *,'      nsound      =  ',nsound
            print *,'      max_nsound  =  ',max_nsound
            print *
            print *,'  Time step (dtl) needs to be smaller '
            print *
            print *,'  -------------------------------- '
            endif
#ifdef MPI
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
            call stopcm1
          endif
        endif

      ENDIF

      end subroutine dtsmall


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcvm(cvm,q)

      use input, only : ib,ie,jb,je,kb,ke,ibm,iem,jbm,jem,kbm,kem,numq,ni,nj,nk,  &
                        eqtset,imoist,iice,nqv,nql1,nql2,nqs1,nqs2,               &
                        timestats,mytime,time_misc
      use constants, only : cv,cvv,cpl,cpi

      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: cvm
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q

      integer :: i,j,k,n

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
    DO k=1,nk

      IF( eqtset.le.1 .or. imoist.eq.0 )THEN

        do j=1,nj
        do i=1,ni
          cvm(i,j,k) = cv
        enddo
        enddo

      ELSE

        do j=1,nj
        do i=1,ni
          cvm(i,j,k) = cv+cvv*q(i,j,k,nqv)
        enddo
        enddo
      IF( nql1.gt.0 )THEN
        do n=nql1,nql2
          do j=1,nj
          do i=1,ni
            cvm(i,j,k)=cvm(i,j,k)+cpl*q(i,j,k,n)
          enddo
          enddo
        enddo
      ENDIF
        IF(iice.eq.1)THEN
          do n=nqs1,nqs2
          do j=1,nj
          do i=1,ni
            cvm(i,j,k)=cvm(i,j,k)+cpi*q(i,j,k,n)
          enddo
          enddo
          enddo
        ENDIF

      ENDIF

    ENDDO

      if(timestats.ge.1) time_misc=time_misc+mytime()

      end subroutine getcvm


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefq(rmax,asq,ruh,rvh,rmh,rho,q3d)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,dx,dy,dz,  &
                        pdscheme,timestats,mytime,time_misc

      implicit none

      real rmax
      double precision :: asq
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,q3d

      integer i,j,k
      double precision :: t1,t2,t3
      double precision :: a1,a2,tem
      double precision, dimension(nj) :: budj
      double precision, dimension(nk) :: budk

!----------------------------------------------------------------------

      tem = dx*dy*dz

      IF(pdscheme.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,t1,t2,t3,a1,a2)
        do j=1,nj
        budj(j)=0.0d0
        do i=1,ni
          t1=0.0d0
          t2=0.0d0
          a1=0.0d0
          a2=0.0d0
          do k=1,nk
            t1=t1+rho(i,j,k)*q3d(i,j,k)
            a1=a1+rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
!!!            q3d(i,j,k)=max(0.0,q3d(i,j,k))
            if(q3d(i,j,k).lt.rmax) q3d(i,j,k)=0.0
            t2=t2+rho(i,j,k)*q3d(i,j,k)
          enddo
          t3=(t1+1.0d-20)/(t2+1.0d-20)
          if(t3.lt.0.0) t3=1.0d0
          do k=1,nk
            q3d(i,j,k)=t3*q3d(i,j,k)
            a2=a2+rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
          enddo
          budj(j)=budj(j)+a2-a1
        enddo
        enddo

        do j=1,nj
          asq=asq+budj(j)*tem
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k,a1,a2)
        do k=1,nk
        budk(k)=0.0d0
        do j=1,nj
        do i=1,ni
          a1=rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
!!!          q3d(i,j,k)=max(0.0,q3d(i,j,k))
          if(q3d(i,j,k).lt.rmax) q3d(i,j,k)=0.0
          a2=rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
          budk(k)=budk(k)+a2-a1
        enddo
        enddo
        enddo

        do k=1,nk
          asq=asq+budk(k)*tem
        enddo

      ENDIF

!----------------------------------------------------------------------

      if(timestats.ge.1) time_misc=time_misc+mytime()

      end subroutine pdefq


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine calcprs(pi0,prs,pp3d)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,  &
                        timestats,mytime,time_prsrho
      use constants, only : p00,cpdrd

      implicit none
 
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0
      real, dimension(ib:ie,jb:je,kb:ke) :: prs,pp3d
 
      integer i,j,k
 
!----------------------------------------------------------------------
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
      enddo
      enddo
      enddo
 
!----------------------------------------------------------------------
 
      if(timestats.ge.1) time_prsrho=time_prsrho+mytime()
 
      end subroutine calcprs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcrho(pi0,th0,rho,prs,pp3d,th3d,q3d)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,               &
                        ibm,iem,jbm,jem,kbm,kem,numq,nqv,imoist,  &
                        timestats,mytime,time_prsrho
      use constants, only : rd,reps

      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,prs,pp3d,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d

      integer i,j,k

!----------------------------------------------------------------------

      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          rho(i,j,k)=prs(i,j,k)                         &
             /( rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
                  *(1.0+max(0.0,q3d(i,j,k,nqv))*reps) )
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------------

      if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

      end subroutine calcrho


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcdbz(rho,qr,qs,qg,dbz)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,  &
                        ptype,dowr,outfile,          &
                        timestats,mytime,time_write
      use goddard_module, only : ROQR,ROQG,ROQS,TNW,TNG,TNSS

      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: rho,qr,qs,qg,dbz

      integer :: i,j,k
      real :: n0r,n0g,n0s,rhor,rhog,rhos,gamma,zer,zeg,zes
      real, parameter :: epp = 1.0e-8

      ! Reference:  Fovell and Ogura, 1988, JAS, pg 3850
      !             (and references therein)

  IF(ptype.eq.2)THEN

    rhor = 1000.0 * ROQR
    rhog = 1000.0 * ROQG
    rhos = 1000.0 * ROQS

    n0r = 1.0e8 * TNW
    n0g = 1.0e8 * TNG
    n0s = 1.0e8 * TNSS

!!!    print *,'  rhor,rhog,rhos = ',rhor,rhog,rhos
!!!    print *,'  n0r,n0g,n0s    = ',n0r,n0g,n0s

!$omp parallel do default(shared)  &
!$omp private(i,j,k,gamma,zer,zeg,zes)
    do k=1,nk
    do j=1,nj
    do i=1,ni

    if(qr(i,j,k).ge.epp)then
      !--- rain ---
      gamma=(3.14159*n0r*rhor/(rho(i,j,k)*qr(i,j,k)))**0.25
      zer=720.0*n0r*(gamma**(-7))
    else
      zer=0.0
    endif

    if(qg(i,j,k).ge.epp)then
      !--- graupel/hail ---
      gamma=(3.14159*n0g*rhog/(rho(i,j,k)*qg(i,j,k)))**0.25
      zeg=720.0*n0g*(gamma**(-7))*((rhog/rhor)**2)*0.224
    else
      zeg=0.0
    endif

    if(qs(i,j,k).ge.epp)then
      !--- snow ---
      gamma=(3.14159*n0s*rhos/(rho(i,j,k)*qs(i,j,k)))**0.25
      zes=720.0*n0s*(gamma**(-7))*((rhos/rhor)**2)*0.224
    else
      zes=0.0
    endif

      !--- dbz ---

    if( (zer+zeg+zes).gt.1.0e-18 )then
      dbz(i,j,k)=10.0*log10((zer+zeg+zes)*1.0e18)
    else
      dbz(i,j,k)=0.0
    endif

    enddo
    enddo
    enddo

  ELSE

    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) ' ptype = ',ptype
    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) ' calcdbz is not valid for this value of ptype'
    if(dowr) write(outfile,*)
    call stopcm1

  ENDIF

      if(timestats.ge.1) time_write=time_write+mytime()

      end subroutine calcdbz


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcuh(uf,vf,zh,zf,ua,va,wa,uh,zeta,dum1,dum2, &
                        zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)

      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ni,nj,nk,  &
                        cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,terrain_flag,  &
                        rdx,rdy,zt

      implicit none

      ! Subroutine to calculate vertically integrated updraft helicity
      ! Reference:  Kain et al, 2008, WAF, p 931

      ! note:  need zh,zf Above Ground Level

      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je) :: uh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zeta,dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: rgzu,rgzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf

      real, parameter :: zz0 = 2000.0     ! bottom of integration layer (m AGL)
      real, parameter :: zzt = 5000.0     ! top of integration layer (m AGL)

      integer :: i,j,k
      real :: r1,r2
      real :: wbar,zbar

  IF(.not.terrain_flag)THEN

    ! Cartesian grid, without terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    DO k=1,nk
    DO j=1,nj+1
    DO i=1,ni+1
      zeta(i,j,k) = (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
                   -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
    ENDDO
    ENDDO
    ENDDO

  ELSE

    ! Cartesian grid, with terrain:

        ! dum1 stores u at w-pts:
        ! dum2 stores v at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+2
          ! lowest model level:
          do i=0,ni+2
            dum1(i,j,1) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
            dum2(i,j,1) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
          enddo

          ! upper-most model level:
          do i=0,ni+2
            dum1(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+2
            dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
            dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
          enddo
          enddo
        enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1)
        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
            zeta(i,j,k)=( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
                         +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                               -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                              )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )        &
                       -( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
                         +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                               -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                              )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )
          enddo
          enddo
        enddo

  ENDIF

!$omp parallel do default(shared)  &
!$omp private(i,j,k,wbar,zbar)
    DO j=1,nj
    DO i=1,ni
      uh(i,j) = 0.0
      DO k=1,nk
        IF( zh(i,j,k).ge.zz0 .and. zh(i,j,k).le.zzt )THEN
          ! note:  only consider cyclonically rotating updrafts
          !        (so, w and zeta must both be positive)
          wbar = max( 0.0 , 0.5*(wa(i,j,k)+wa(i,j,k+1)) )
          zbar = max( 0.0 , 0.25*(zeta(i,j,k)+zeta(i+1,j,k)   &
                                 +zeta(i,j+1,k)+zeta(i+1,j+1,k)) )
          uh(i,j) = uh(i,j) + (min(zf(i,j,k+1),zzt)-max(zf(i,j,k),zz0))*wbar*zbar
        ENDIF
      ENDDO
    ENDDO
    ENDDO

      end subroutine calcuh


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcvort(xh,xf,uf,vf,zh,mh,zf,mf,                                         &
                          zs,gz,gzu,gzv,rgz,rgzu,rgzv,gxu,gyv,rds,sigma,rdsf,sigmaf,       &
                          ua,va,wa,xvort,yvort,zvort,tem ,dum1,dum2,pv  ,th  ,th0,tha,rr,  &
                          ust,znt,u1,v1,s1,do_avg)

      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ktb,kte,ni,nj,nk,  &
                        cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,terrain_flag,          &
                        dgs1,dgs2,dgs3,dgt1,dgt2,dgt3,imove,umove,vmove,     &
                        rdx,rdy,rdz,zt,bbc,tbc,axisymm,output_pv,myid
      use constants

      implicit none

      ! Subroutine to calculate 3 components of vorticity
      ! at scalar points.

      ! cm1r19.6:  ua and va are now ground-relative winds when imove=1

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,gzu,gzv,rgz,rgzu,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu,gyv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: xvort,yvort,zvort,tem,dum1,dum2,pv,th
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,tha,rr
      real, intent(in), dimension(ib:ie,jb:je) :: ust,znt,u1,v1,s1
      logical, intent(in), optional :: do_avg

      integer :: i,j,k
      real :: r1,r2,ustavg,uavg,vavg,savg,zavg
      logical :: doavg

      doavg = .true.
      if( present(do_avg) )then
        if( .not.do_avg ) doavg = .false.
      endif

      if(myid.eq.0) print *,'  doavg = ',doavg

!-----------------------------------------------------------------------

      IF( output_pv.eq.1 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pv(i,j,k)=0.0
        enddo
        enddo
        enddo
      ENDIF

    IF( terrain_flag .or. output_pv.eq.1 )THEN
      ! dum1 stores w at scalar pts
      ! dum2 stores theta at w pts
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do j=0,nj+1
        do k=1,nk
        do i=0,ni+1
          dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
          th(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        ! lowest model level:
        do i=0,ni+1
          dum2(i,j,1) = cgs1*th(i,j,1)+cgs2*th(i,j,2)+cgs3*th(i,j,3)
        enddo
        ! upper-most model level:
        do i=0,ni+1
          dum2(i,j,nk+1) = cgt1*th(i,j,nk)+cgt2*th(i,j,nk-1)+cgt3*th(i,j,nk-2)
        enddo
        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+1
          dum2(i,j,k) = r1*th(i,j,k-1)+r2*th(i,j,k)
        enddo
        enddo
      enddo
    ENDIF

!-----------------------------------------------------------------------
! x-vort:

  tem=0.0
  if(axisymm.eq.0)then
    IF(.not.terrain_flag)THEN
      !cccccccccccccccccccccccccccccccccccccc
      ! Cartesian grid, without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj+1
      do i=1,ni
        tem(i,j,k) = (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                    -(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
      enddo
      enddo
      enddo
      IF( bbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,1)=tem(i,j,2)
        enddo
        enddo
      ELSEIF( bbc.eq.2 )THEN
      if( imove.eq.0 )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,1)=-2.0*va(i,j,1)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
        enddo
        enddo
      else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,1)=-2.0*(va(i,j,1)+vmove)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
        enddo
        enddo
      endif
      ENDIF
      IF( tbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,nk+1)=tem(i,j,nk)
        enddo
        enddo
      ELSEIF( tbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,nk+1)=2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
        enddo
        enddo
      ENDIF
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          xvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i,j+1,k)+tem(i,j,k+1)+tem(i,j+1,k+1))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj+1
        do i=1,ni
          xvort(i,j,k) = tem(i,j,k)
        enddo
        enddo
        enddo
      endif
      !cccccccccccccccccccccccccccccccccccccc
    getpv11: IF( output_pv.eq.1 )THEN
      ! here, zvort array stores d(th)/dx
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        zvort(i,j,k) = (th(i,j,k)-th(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo
      ! pv1:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pv(i,j,k)=pv(i,j,k)+xvort(i,j,k)*0.5*( zvort(i,j,k)+zvort(i+1,j,k) )
      enddo
      enddo
      enddo
    ENDIF  getpv11
      !cccccccccccccccccccccccccccccccccccccc
    ELSE
      !cccccccccccccccccccccccccccccccccccccc
      ! Cartesian grid, with terrain:
      !   (dum1 stores w at scalar-pts:)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj+1
      do i=1,ni
        tem(i,j,k)=(-(va(i,j,k)-va(i,j,k-1))*rds(k)                                  &
                    +(wa(i,j,k)*rgz(i,j)-wa(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)           &
                    +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))    &
                                 -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )  &
                               *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)                      &
                   )*0.5*( gz(i,j)+gz(i,j-1) )
      enddo
      enddo
      enddo
      IF( bbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,1)=tem(i,j,2)
        enddo
        enddo
      ELSEIF( bbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,1)=-2.0*va(i,j,1)*rds(2)*0.5*( gz(i,j-1)+gz(i,j) )
        enddo
        enddo
      ENDIF
      IF( tbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,nk+1)=tem(i,j,nk)
        enddo
        enddo
      ELSEIF( tbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni
          tem(i,j,nk+1)=2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
        enddo
        enddo
      ENDIF
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          xvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i,j+1,k)+tem(i,j,k+1)+tem(i,j+1,k+1))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj+1
        do i=1,ni
          xvort(i,j,k) = tem(i,j,k)
        enddo
        enddo
        enddo
      endif
    getpv1: IF( output_pv.eq.1 )THEN
      !cccccccccccccccccccccccccccccccccccccc
      ! here, zvort array stores d(th)/dx
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        zvort(i,j,k) = gzu(i,j)*(th(i,j,k)*rgz(i,j)-th(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)  &
               +0.5*( gxu(i,j,k+1)*(dum2(i,j,k+1)+dum2(i-1,j,k+1))                     &
                     -gxu(i,j,k  )*(dum2(i,j,k  )+dum2(i-1,j,k  )) )*rdsf(k)
      enddo
      enddo
      enddo
      ! pv1:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pv(i,j,k)=pv(i,j,k)+xvort(i,j,k)*0.5*( zvort(i,j,k)+zvort(i+1,j,k) )
      enddo
      enddo
      enddo
      !cccccccccccccccccccccccccccccccccccccc
    ENDIF  getpv1
    ENDIF
  else
      !cccccccccccccccccccccccccccccccccccccc
      ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tem(i,j,k) = -(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo
      enddo
      IF( bbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          tem(i,j,1)=tem(i,j,2)
        enddo
        enddo
      ELSEIF( bbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          tem(i,j,1)=-2.0*va(i,j,1)*rdz*mf(1,1,1)
        enddo
        enddo
      ENDIF
      IF( tbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          tem(i,j,nk+1)=tem(i,j,nk)
        enddo
        enddo
      ELSEIF( tbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          tem(i,j,nk+1)=2.0*va(i,j,nk)*rdz*mf(1,1,nk+1)
        enddo
        enddo
      ENDIF
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          xvort(i,j,k) = 0.5*(tem(i,j,k)+tem(i,j,k+1))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj+1
        do i=1,ni
          xvort(i,j,k) = tem(i,1,k)
        enddo
        enddo
        enddo
      endif
      !cccccccccccccccccccccccccccccccccccccc
  endif

!-----------------------------------------------------------------------
! y-vort:

    tem=0.0
    IF(.not.terrain_flag)THEN
      !cccccccccccccccccccccccccccccccccccccc
      ! Cartesian grid, without terrain:
      ! and axisymmetric grid:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni+1
        tem(i,j,k) = (ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))   &
                    -(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo
      IF( bbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,1)=tem(i,j,2)
        enddo
        enddo
      ELSEIF( bbc.eq.2 )THEN
      if( imove.eq.0 )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,1)=2.0*ua(i,j,1)*rdz*0.5*(mf(i-1,j,1)+mf(i,j,1))
        enddo
        enddo
      else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,1)=2.0*(ua(i,j,1)+umove)*rdz*0.5*(mf(i-1,j,1)+mf(i,j,1))
        enddo
        enddo
      endif
      ENDIF
      IF( tbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,nk+1)=tem(i,j,nk)
        enddo
        enddo
      ELSEIF( tbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
        enddo
        enddo
      ENDIF
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          yvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j,k+1)+tem(i+1,j,k+1))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni+1
          yvort(i,j,k) = tem(i,j,k)
        enddo
        enddo
        enddo
      endif
      !cccccccccccccccccccccccccccccccccccccc
    getpv12:  IF( output_pv.eq.1 )THEN
      ! here, zvort array stores d(th)/dy
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        zvort(i,j,k) = (th(i,j,k)-th(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
      ! pv1:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pv(i,j,k)=pv(i,j,k)+yvort(i,j,k)*0.5*( zvort(i,j,k)+zvort(i,j+1,k) )
      enddo
      enddo
      enddo
    ENDIF  getpv12
      !cccccccccccccccccccccccccccccccccccccc
    ELSE
      !cccccccccccccccccccccccccccccccccccccc
      ! Cartesian grid, with terrain:
      !   (dum1 stores w at scalar-pts:)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni+1
        tem(i,j,k)=( (ua(i,j,k)-ua(i,j,k-1))*rds(k)                                  &
                    -(wa(i,j,k)*rgz(i,j)-wa(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)           &
                    -0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))    &
                                 -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )  &
                               *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)                      &
                   )*0.5*( gz(i,j)+gz(i-1,j) )
      enddo
      enddo
      enddo
      IF( bbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,1)=tem(i,j,2)
        enddo
        enddo
      ELSEIF( bbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,1)=2.0*ua(i,j,1)*rds(2)*0.5*( gz(i-1,j)+gz(i,j) )
        enddo
        enddo
      ENDIF
      IF( tbc.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,nk+1)=tem(i,j,nk)
        enddo
        enddo
      ELSEIF( tbc.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni+1
          tem(i,j,nk+1)=-2.0*ua(i,j,nk)*rds(nk)*0.5*( gz(i-1,j)+gz(i,j) )
        enddo
        enddo
      ENDIF
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          yvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j,k+1)+tem(i+1,j,k+1))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni+1
          yvort(i,j,k) = tem(i,j,k)
        enddo
        enddo
        enddo
      endif
    getpv2:  IF( output_pv.eq.1 )THEN
      !cccccccccccccccccccccccccccccccccccccc
      ! here, zvort array stores d(th)/dy
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        zvort(i,j,k) = gzv(i,j)*(th(i,j,k)*rgz(i,j)-th(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)  &
               +0.5*( gyv(i,j,k+1)*(dum2(i,j,k+1)+dum2(i,j-1,k+1))                     &
                     -gyv(i,j,k  )*(dum2(i,j,k  )+dum2(i,j-1,k  )) )*rdsf(k)
      enddo
      enddo
      enddo
      ! pv1:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pv(i,j,k)=pv(i,j,k)+yvort(i,j,k)*0.5*( zvort(i,j,k)+zvort(i,j+1,k) )
      enddo
      enddo
      enddo
    ENDIF  getpv2
      !cccccccccccccccccccccccccccccccccccccc
    ENDIF

!-----------------------------------------------------------------------
! z-vort:

    tem=0.0
    if(axisymm.eq.0)then
      IF(.not.terrain_flag)THEN
        !cccccccccccccccccccccccccccccccccccccc
        ! Cartesian grid, without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            tem(i,j,k) = (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
                        -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
          enddo
          enddo
        enddo
        if( doavg )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            zvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j+1,k)+tem(i+1,j+1,k))
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            zvort(i,j,k) = tem(i,j,k)
          enddo
          enddo
          enddo
        endif
        !cccccccccccccccccccccccccccccccccccccc
      ELSE
        !cccccccccccccccccccccccccccccccccccccc
        ! Cartesian grid, with terrain:
        ! dum1 stores u at w-pts:
        ! dum2 stores v at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+2
          ! lowest model level:
          do i=0,ni+2
            dum1(i,j,1) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
            dum2(i,j,1) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
          enddo
          ! upper-most model level:
          do i=0,ni+2
            dum1(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
          enddo
          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+2
            dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
            dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
          enddo
          enddo
        enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1)
        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
            tem(i,j,k)=r1*( (va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
                           +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))   &
                                 -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))   &
                                )*rdsf(k)*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )        &
                      -r1*( (ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
                           +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))   &
                                 -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))   &
                                )*rdsf(k)*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )
          enddo
          enddo
        enddo
        if( doavg )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            zvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j+1,k)+tem(i+1,j+1,k))
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            zvort(i,j,k) = tem(i,j,k)
          enddo
          enddo
          enddo
        endif
        !cccccccccccccccccccccccccccccccccccccc
      ENDIF
    else
      !cccccccccccccccccccccccccccccccccccccc
      !$omp parallel do default(shared)  &
      !$omp private(j,k)
      do k=1,nk
      do j=1,nj
        tem(1,j,k) = 0.0
      enddo
      enddo
      ! Axisymmetric grid:
      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        tem(i,j,k) = (va(i,j,k)*xh(i)-va(i-1,j,k)*xh(i-1))*rdx*uf(i)/xf(i)
      enddo
      enddo
      enddo
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          zvort(i,j,k) = 0.5*(tem(i,j,k)+tem(i+1,j,k))
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          zvort(i,j,k) = tem(i,1,k)
        enddo
        enddo
        enddo
      endif
      !cccccccccccccccccccccccccccccccccccccc
    endif


    getpv3:  IF( output_pv.eq.1 )THEN
      ! now, dum1 stores dt/dz:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do j=1,nj
        do k=2,nk
        do i=1,ni
          dum1(i,j,k) = (th(i,j,k)-th(i,j,k-1))*rdz*mf(i,j,k)
        enddo
        enddo
        do i=1,ni
          dum1(i,j,1) = (dgs3*th(i,j,3)+dgs2*th(i,j,2)+dgs1*th(i,j,1))*rdz*mh(i,j,1)
        enddo
        do i=1,ni
          dum1(i,j,nk+1) = (dgt3*th(i,j,nk-2)+dgt2*th(i,j,nk-1)+dgt1*th(i,j,nk))*rdz*mh(i,j,nk)
        enddo
        ! pv:
        do k=1,nk
        do i=1,ni
          pv(i,j,k)=pv(i,j,k)+zvort(i,j,k)*0.5*(dum1(i,j,k)+dum1(i,j,k+1))
          pv(i,j,k)=pv(i,j,k)*rr(i,j,k)
        enddo
        enddo
      enddo
    ENDIF  getpv3

!-----------------------------------------------------------------------

    IF( bbc.eq.3 )THEN
      ! cm1r18:  use log-layer equation below
      if( doavg )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          xvort(i,j,1) = -(ust(i,j)/(karman*zh(i,j,1)))*(v1(i,j)/max(s1(i,j),0.01))
          yvort(i,j,1) =  (ust(i,j)/(karman*zh(i,j,1)))*(u1(i,j)/max(s1(i,j),0.01))
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,ustavg,uavg,vavg,savg,zavg)
        do j=1,nj
        do i=1,ni
          ustavg = 0.5*(ust(i,j-1)+ust(i,j))
          vavg = 0.5*(v1(i,j-1)+v1(i,j))
          savg = 0.5*(s1(i,j-1)+s1(i,j))
          zavg = 0.5*(zh(i,j-1,1)+zh(i,j,1))
          xvort(i,j,1) = -(ustavg/(karman*zavg))*(vavg/max(savg,0.01))
          ustavg = 0.5*(ust(i-1,j)+ust(i,j))
          uavg = 0.5*(v1(i-1,j)+v1(i,j))
          savg = 0.5*(s1(i-1,j)+s1(i,j))
          zavg = 0.5*(zh(i-1,j,1)+zh(i,j,1))
          yvort(i,j,1) =  (ustavg/(karman*zavg))*(uavg/max(savg,0.01))
        enddo
        enddo
      endif
    ENDIF

!-----------------------------------------------------------------------

      end subroutine calcvort


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calccpch(zh,zf,th0,qv0,cpc,cph,tha,qa)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,        &
                        ibm,iem,jbm,jem,kbm,kem,numq,nqv,  &
                        imoist,iice,nql1,nql2,nqs1,nqs2
      use constants, only : g,repsm1

      implicit none

      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: th0,qv0
      real, intent(inout), dimension(ib:ie,jb:je) :: cpc,cph
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: tha
      real, intent(in),    dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa

      integer :: i,j,k,n
      real :: ql
      real, dimension(nk) :: bb

      ! defines top of cold pool / location to stop calculation of C
      real, parameter :: bcrit = -0.01

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,ql,bb)
    DO j=1,nj
    DO i=1,ni
      cpc(i,j) = 0.0
      cph(i,j) = 0.0
      bb = 0.0
      do k=1,nk
        bb(k) = g*tha(i,j,k)/th0(i,j,k)
      enddo
      if(imoist.eq.1)then
        do k=1,nk
          ql = 0.0
        IF( nql1.gt.0 )THEN
          do n=nql1,nql2
            ql=ql+qa(i,j,k,n)
          enddo
        ENDIF
          if(iice.eq.1)then
            do n=nqs1,nqs2
              ql=ql+qa(i,j,k,n)
            enddo
          endif
          bb(k) = bb(k) + g*( repsm1*(qa(i,j,k,nqv)-qv0(i,j,k)) - ql )
        enddo
      endif
    ! only calculate cpc/cph if surface B is less than bcrit
    IF( bb(1).lt.bcrit .and. tha(i,j,1).le.-1.0 )THEN
      cpc(i,j) = - 2.0*bb(1)*(zf(i,j,2)-zf(i,j,1))
      k = 2
      do while( bb(k).lt.bcrit .and. k.lt.nk )
        if( cpc(i,j).lt.0.0 ) cpc(i,j) = 0.0
        cpc(i,j) = cpc(i,j) - 2.0*bb(k)*(zf(i,j,k+1)-zf(i,j,k))
        k = k + 1
      enddo
      if( cpc(i,j).gt.0.0 )then
        cpc(i,j) = sqrt(cpc(i,j))
        if(k.eq.nk)then
          cph(i,j) = zf(i,j,nk+1)
        else
          cph(i,j) = zh(i,j,k-1) + (zh(i,j,k)-zh(i,j,k-1))*(bcrit-bb(k-1))   &
                                                          /(bb(k)-bb(k-1))
        endif
        ! account for terrain:
        cph(i,j) = cph(i,j) - zf(i,j,1)
      endif
    ENDIF
    ENDDO
    ENDDO

      end subroutine calccpch


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calccref(cref,dbz)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,  &
                        timestats,mytime,time_write

      implicit none

      real, intent(inout), dimension(ib:ie,jb:je) :: cref
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: dbz

      integer :: i,j,k

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        cref(i,j) = -1000.0
      enddo
      enddo

      do k=1,nk

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        cref(i,j)=max(cref(i,j),dbz(i,j,k))
      enddo
      enddo

      enddo


      if(timestats.ge.1) time_write=time_write+mytime()

      end subroutine calccref


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calccflquick(dt,uh,vh,mh,ua,va,wa,reqc)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,nx,ny,                 &
                        rdx,rdy,rdz,axisymm,cflmax,myid,                  &
                        timestats,mytime,time_cflq

#ifdef MPI
      use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      integer, intent(inout) :: reqc
 
      integer :: i,j,k,ierr
      real :: dtdx,dtdy,dtdz,fmax,wsp
      real, dimension(nk) :: cfl
 
      dtdx=0.5*dt*rdx
      dtdy=0.5*dt*rdy
      dtdz=0.5*dt*rdz

      !$omp parallel do default(shared)  &
      !$omp private(k)
      do k=1,nk
        cfl(k) = 0.0
      enddo

      if(nx.gt.1.and.ny.gt.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,wsp)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          wsp = sqrt( ( ((ua(i,j,k)+ua(i+1,j,k))*dtdx*uh(i))**2     &
                       +((va(i,j,k)+va(i,j+1,k))*dtdy*vh(j))**2 )   &
                       +((wa(i,j,k)+wa(i,j,k+1))*dtdz*mh(i,j,k))**2 )
          cfl(k) = max( cfl(k) , wsp )
!!!          if( wsp.ge.1.30 ) print *,'  cfl,myid,i,j,k = ',wsp,myid,i,j,k,ua(i,j,k),ua(i+1,j,k),va(i,j,k),va(i,j+1,k),wa(i,j,k),wa(i,j,k+1)
        enddo
        enddo
        enddo
      elseif(nx.gt.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,wsp)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          wsp = sqrt( ((ua(i,j,k)+ua(i+1,j,k))*dtdx*uh(i))**2     &
                     +((wa(i,j,k)+wa(i,j,k+1))*dtdz*mh(i,j,k))**2 )
          cfl(k) = max( cfl(k) , wsp )
!!!          if( wsp.ge.1.30 ) print *,'  cfl,myid,i,j,k = ',wsp,myid,i,j,k
        enddo
        enddo
        enddo
      elseif(axisymm.eq.0.and.ny.gt.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,wsp)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          wsp = sqrt( ((va(i,j,k)+va(i,j+1,k))*dtdy*vh(j))**2     &
                     +((wa(i,j,k)+wa(i,j,k+1))*dtdz*mh(i,j,k))**2 )
          cfl(k) = max( cfl(k) , wsp )
!!!          if( wsp.ge.1.30 ) print *,'  cfl,myid,i,j,k = ',wsp,myid,i,j,k
        enddo
        enddo
        enddo
      endif

      fmax=-99999999.
      do k=1,nk
        fmax = max( fmax , cfl(k) )
      enddo

      cflmax = fmax

#ifdef MPI
      call MPI_IALLREDUCE(MPI_IN_PLACE,cflmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,reqc,ierr)
#endif

      if(timestats.ge.1) time_cflq=time_cflq+mytime()
 
      end subroutine calccflquick


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcksquick(dt,uh,vh,mf,kmh,kmv,khh,khv,reqk,cm0,xh,yh,zh)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,nx,ny,nz,      &
                        rdx,rdy,rdz,doimpl,cm1setup,ksmax,        &
                        ibc,iec,jbc,jec,kbc,kec,axisymm,          &
                        timestats,mytime,time_cflq
      use constants, only : cmemin
#ifdef MPI
      use mpi
#endif

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      integer, intent(inout) :: reqk
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
 
      integer :: i,j,k,ierr
      real :: dtdx,dtdy,dtdz,fmax,tem1,tem2,tem3,tem4
      real, dimension(2:nk) :: ks
 
      dtdx=dt*rdx*rdx
      dtdy=dt*rdy*rdy
      if( doimpl.eq.0 )then
        dtdz=dt*rdz*rdz
      else
        dtdz=0.0
      endif

      ks = 0.0

      if(nx.gt.1.and.ny.gt.1)then
        ! 3d:
      if( cm1setup.eq.4 )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,tem1,tem2,tem3,tem4)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          ! LES within MM 
        if( cm0(i  ,j  ).gt.cmemin .or.  &
            cm0(i+1,j  ).gt.cmemin .or.  &
            cm0(i-1,j  ).gt.cmemin .or.  &
            cm0(i  ,j+1).gt.cmemin .or.  &
            cm0(i  ,j-1).gt.cmemin )then
          tem1 = sqrt( (kmh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(kmh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem2 = sqrt( (khh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(khh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem3 = sqrt( (kmv(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(kmv(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem4 = max(kmv(i,j,k),khv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k)
          ks(k) = max( ks(k) , tem1 , tem2 , tem3 , tem4 )
        else
          tem1 = sqrt( (kmh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(kmh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem2 = sqrt( (khh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(khh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem4 = max(kmv(i,j,k),khv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k)
          ks(k) = max( ks(k) , tem1 , tem2 , tem4 )
        endif
        enddo
        enddo
        enddo
      else
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,tem1,tem2,tem3)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          ! check entire domain
          tem1 = sqrt( (kmh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(kmh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem2 = sqrt( (khh(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(khh(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem3 = sqrt( (kmv(i,j,k)*dtdx*uh(i)*uh(i))**2 &
                      +(kmv(i,j,k)*dtdy*vh(j)*vh(j))**2 )
          tem4 = max(kmv(i,j,k),khv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k)
          ks(k) = max( ks(k) , tem1 , tem2 , tem3 , tem4 )
        enddo
        enddo
        enddo
      endif
      elseif(nx.gt.1)then
        ! 2d (including axisymm):
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,tem1,tem2,tem3)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tem1 = kmh(i,j,k)*dtdx*uh(i)*uh(i)
          tem2 = khh(i,j,k)*dtdx*uh(i)*uh(i)
          tem3 = kmv(i,j,k)*dtdx*uh(i)*uh(i)
          tem4 = max(kmv(i,j,k),khv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k)
          ks(k) = max( ks(k) , tem1 , tem2 , tem3 , tem4 )
        enddo
        enddo
        enddo
      elseif(axisymm.eq.0.and.ny.gt.1)then
        stop 1112
      endif

      fmax=-99999999.
      do k=2,nk
        fmax = max( fmax , ks(k) )
      enddo

      ksmax = fmax

#ifdef MPI
      call MPI_IALLREDUCE(MPI_IN_PLACE,ksmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,reqk,ierr)
#endif

      if(timestats.ge.1) time_cflq=time_cflq+mytime()
      end subroutine calcksquick


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine check_for_stop()

      use input, only : stopit,stopcondition,cflmax,smeps,myid
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

      implicit none

      if( cflmax.ge.1.50 )then
        ! this value (1.5) is arbitrary; could potentialy be higher, 
        ! maybe should be lower
        stopit=.true.
        stopcondition = 1
      elseif( cflmax.lt.-smeps )then
        stopit=.true.
        stopcondition = 2
      elseif( ieee_is_nan(cflmax) )then
        stopit=.true.
        stopcondition = 3
      endif

      if(stopit)then
        if(myid.eq.0) print *,'  cflmax = ',cflmax
      endif

      end subroutine check_for_stop


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine zinterp(sigma,zs,zh,dum1,dum2)

      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk
      use constants, only : grads_undef

      implicit none

      real, dimension(kb:ke) :: sigma
      real, dimension(ib:ie,jb:je) :: zs
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,dum1,dum2

      integer i,j,k,kk,kup,kdn
      real, dimension(nk) :: zref

      do k=1,nk
!!!        zref(k)=(k*dz-0.5*dz)
        zref(k)=sigma(k)
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum2(i,j,k)=dum1(i,j,k)
      enddo
      enddo
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        if( (zref(k).lt.zh(i,j,1)).or.(zref(k).gt.zh(i,j,nk)) )then
          if( zref(k).gt.0.5*zh(i,j,1) .and. zref(k).gt.zs(i,j) )then
            ! 2nd-order extrapolation:
            dum1(i,j,k)=dum2(i,j,1)-(zh(i,j,1)-zref(k))                             &
                                   *(-3.0*dum2(i,j,1)+4.0*dum2(i,j,2)-dum2(i,j,3))  &
                                   *0.25/(zh(i,j,1)-zs(i,j))
          else
            dum1(i,j,k)=grads_undef
          endif
        elseif(zs(i,j).lt.0.1 .or. zref(k).eq.zh(i,j,1))then
          dum1(i,j,k)=dum2(i,j,k)
        else
          kup=0
          kdn=0
          do kk=1,nk
            if(zref(k).gt.zh(i,j,kk)) kdn=kk
          enddo
          kup=kdn+1
          if(kup.le.0.or.kdn.le.0.or.kup.ge.nk+1.or.kdn.ge.nk+1)then
            print *,kdn,kup
            print *,zs(i,j),zh(i,j,kdn),zref(k),zh(i,j,kup)
            print *,i,j,k
            call stopcm1
          endif
          dum1(i,j,k)=dum2(i,j,kdn)+(dum2(i,j,kup)-dum2(i,j,kdn))   &
                                   *(  zref(k  )  -zh(i,j,kdn))     &
                                   /(  zh(i,j,kup)-zh(i,j,kdn))
        endif
      enddo
      enddo
      enddo

      end subroutine zinterp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine set_time_to_zero
      use input
      implicit none

      time_sound=0.0
      time_poiss=0.0
      time_advs=0.0
      time_advu=0.0
      time_advv=0.0
      time_advw=0.0
      time_buoyan=0.0
      time_turb=0.0
      time_diffu=0.0
      time_microphy=0.0
      time_dbz=0.0
      time_stat=0.0
      time_cflq=0.0
      time_bc=0.0
      time_misc=0.0
      time_integ=0.0
      time_rdamp=0.0
      time_divx=0.0
      time_write=0.0
      time_restart=0.0
      time_ttend=0.0
      time_cor=0.0
      time_fall=0.0
      time_satadj=0.0
      time_sfcphys=0.0
      time_parcels=0.0
      time_hailtrajs=0.0
      time_rad=0.0
      time_pbl=0.0
      time_swath=0.0
      time_pdef=0.0
      time_prsrho=0.0
      time_diag=0.0
      time_azimavg=0.0
      time_hifrq=0.0
      time_ercyl=0.0
      time_tavg=0.0
      time_cuprm=0.0
      time_param=0.0
#ifdef MPI
      time_mpu1=0.0
      time_mpv1=0.0
      time_mpw1=0.0
      time_mpp1=0.0
      time_mpu2=0.0
      time_mpv2=0.0
      time_mpw2=0.0
      time_mpp2=0.0
      time_mps1=0.0
      time_mps3=0.0
      time_mpq1=0.0
      time_mptk1=0.0
      time_mptk2=0.0
      time_mps2=0.0
      time_mps4=0.0
      time_mpq2=0.0
      time_mpb=0.0
#endif

      end subroutine set_time_to_zero


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE misclibs
