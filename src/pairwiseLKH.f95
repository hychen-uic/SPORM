!######################################################################

    module pwLKH
    use, intrinsic :: iso_c_binding
    implicit none
    public :: pwMLECOL, pwMLErow

    contains

!######################################################################
!-------------------------------------------
!  pairwise likelihood approach
!  parameter gam is vectorized by column.
!-------------------------------------------

subroutine pwMLECOL(y,x,n,p,q,theta,estv,niter,eps,converge) bind(C, name = "pwmlecol_")

  implicit none
  !integer, parameter:: dp = selected_real_kind(15, 307)

  integer (C_INT) n,p,q,niter,converge
  real (C_DOUBLE)  y(n,p),x(n,q)
  real (C_DOUBLE)  theta(p*q),estv(p*q,p*q) ! columnwise vectorization

  integer (C_INT) i,j,k,kk,l,ll,ki,li,irep
  real (C_DOUBLE)  temp(q,p),probtemp
  real (C_DOUBLE)  der(p*q),der2(p*q,p*q),delta(p*q),eps

  real (C_DOUBLE)  tder(n,p*q)

  converge=0
  do irep=1,niter
    !write(*,*)irep,niter
    der=0
    der2=0
    do i=1,n
      do j=i+1,n

        probtemp=0
        do k=1,q
          do kk=1,p
            temp(k,kk)=(y(i,kk)-y(j,kk))*(x(j,k)-x(i,k))
            probtemp=probtemp+theta((k-1)*p+kk)*temp(k,kk)
          enddo
        enddo
        if(probtemp<=0) then
          probtemp=exp(probtemp)/(1+exp(probtemp))
        else
          probtemp=1.0/(1+exp(-probtemp))
        endif

        do k=1,q
          do kk=1,p
            ki=(k-1)*p+kk
            der(ki)=der(ki)+temp(k,kk)*probtemp
            do l=1,k
              do ll=1,p
                li=(l-1)*p+ll
                if(li<=ki) then
                  der2(ki,li)=der2(ki,li)+temp(k,kk)*temp(l,ll)*probtemp*(1.0-probtemp)
                endif
              enddo
            enddo

          enddo
        enddo

      enddo
    enddo

    do k=1,p*q
      do kk=1,k
        der2(kk,k)=der2(k,kk)
      enddo
    enddo

    Call PDmatinv(der2,p*q)

    delta=matmul(der2,der)
    if(sum(abs(delta))<eps) then
      converge=1
      ! compute variance of the estimator using sandwich estimator

      do i=1,n

        tder(i,:)=0
        do j=1,n
          probtemp=0
          do k=1,q
            do kk=1,p
              temp(k,kk)=(y(i,kk)-y(j,kk))*(x(i,k)-x(j,k))
              probtemp=probtemp+theta((k-1)*p+kk)*temp(k,kk)
            enddo
          enddo
          if(probtemp<=0) then
            probtemp=exp(probtemp)/(1+exp(probtemp))
          else
            probtemp=1.0/(1+exp(-probtemp))
          endif

          do k=1,q
            do kk=1,p
              ki=(k-1)*p+kk
              tder(i,ki)=tder(i,ki)+temp(k,kk)*(1.0-probtemp)
            enddo
          enddo

        enddo
        tder(i,:)=tder(i,:)/(n-1.0)

      enddo

      do k=1,p*q
        do kk=1,k
          estv(k,kk)=0
          do i=1,n
            estv(k,kk)=estv(k,kk)+tder(i,k)*tder(i,kk)
          enddo
          estv(k,kk)=estv(k,kk)*4/(n-1.0)

          estv(kk,k)=estv(k,kk)
        enddo
      enddo

      der2=der2*n*(n-1)/2.0
      estv=matmul(der2,matmul(estv,der2)) ! sandwich estimate of variance

      estv=estv/n ! variance for the estimate
      exit
    else
      theta=theta-delta
    endif

  enddo

end subroutine pwMLECOL
!-------------------------------------------
!  pairwise likelihood approach
!  parameter gam is vectorized by row.
!-------------------------------------------
subroutine pwMLErow(y,x,n,p,q,theta,estv,niter,eps,converge) bind(C, name = "pwmlerow_")

  implicit none
  ! integer, parameter:: dp = selected_real_kind(15, 307)

  integer (C_INT) converge !=1 if converges, =0 otherwise.
  integer (C_INT) n,p,q,niter
  real (C_DOUBLE)  y(n,p),x(n,q)
  real (C_DOUBLE)  theta(p*q),estv(p*q,p*q) ! columnwise vectorization

  integer (C_INT) i,j,k,kk,l,ll,ki,li,irep
  real (C_DOUBLE)  temp(p,q),probtemp
  real (C_DOUBLE)  der(p*q),der2(p*q,p*q),delta(p*q),eps

  real (C_DOUBLE)  tder(n,p*q)

  converge=0
  do irep=1,niter
    !write(*,*)irep,niter
    der=0
    der2=0
    do i=1,n
      do j=i+1,n

        probtemp=0
        do k=1,p
          do kk=1,q
            temp(k,kk)=(y(i,k)-y(j,k))*(x(j,kk)-x(i,kk))
            probtemp=probtemp+theta((k-1)*q+kk)*temp(k,kk)
          enddo
        enddo
        if(probtemp<=0) then
          probtemp=exp(probtemp)/(1+exp(probtemp))
        else
          probtemp=1.0/(1+exp(-probtemp))
        endif

        do k=1,p
          do kk=1,q
            ki=(k-1)*q+kk
            der(ki)=der(ki)+temp(k,kk)*probtemp
            do l=1,k
              do ll=1,q
                li=(l-1)*q+ll
                if(li<=ki) then
                  der2(ki,li)=der2(ki,li)+temp(k,kk)*temp(l,ll)*probtemp*(1-probtemp)
                endif
              enddo
            enddo

          enddo
        enddo

      enddo
    enddo

    do k=1,p*q
      do kk=1,k
        der2(kk,k)=der2(k,kk)
      enddo
    enddo

    Call PDmatinv(der2,p*q)

    delta=matmul(der2,der)
    if(sum(abs(delta))<eps) then
     converge=1
      ! compute variance of the estimator using sandwich estimator

      do i=1,n

        tder(i,:)=0
        do j=1,n
          probtemp=0
          do k=1,p
            do kk=1,q
              temp(k,kk)=(y(i,kk)-y(j,kk))*(x(i,k)-x(j,k))
              probtemp=probtemp+theta((k-1)*q+kk)*temp(k,kk)
            enddo
          enddo
          if(probtemp<=0) then
            probtemp=exp(probtemp)/(1+exp(probtemp))
          else
            probtemp=1.0/(1+exp(-probtemp))
          endif

          do k=1,p
            do kk=1,q
              ki=(k-1)*q+kk
              tder(i,ki)=tder(i,ki)+temp(k,kk)*(1.0-probtemp)
            enddo
          enddo

        enddo
        tder(i,:)=tder(i,:)/(n-1.0)

      enddo

      do k=1,p*q
        !tempm(k)=sum(tder(:,k))/n
 
        do kk=1,k
          estv(k,kk)=0
          do i=1,n
            estv(k,kk)=estv(k,kk)+tder(i,k)*tder(i,kk)
          enddo
          estv(k,kk)=estv(k,kk)*4/(n-1.0)

          estv(kk,k)=estv(k,kk)
        enddo
      enddo

      der2=der2*n*(n-1)/2.0
      estv=matmul(der2,matmul(estv,der2)) ! sandwich estimate of variance

      estv=estv/n ! variance for the estimate
      exit
    else
      theta=theta-delta
    endif

  enddo

end subroutine pwMLErow

!--------------------------------------------
! Simple inverse by Cholestry decomposition
!---------------------------------------------
subroutine PDmatinv(a,n)

      implicit none
      integer, parameter:: dp = selected_real_kind(15, 307)

      integer n
      real(kind=dp)  a(n,n)

      integer i,j,k
      real(kind=dp)  tot,b(n)

      b(1)=sqrt(a(1,1))
      do j=2,n
        a(j,1)=a(j,1)/b(1)
      enddo

      do i=2,n
         do j=i,n
            tot=a(i,j)
            do k=i-1,1,-1
               tot=tot-a(i,k)*a(j,k)
               if(i.eq.j) then
                  if(tot<=0) then
                     write(*,*)j,j,'  choldc failed'
                  endif
                  b(i)=sqrt(tot)
               else
                  a(j,i)=tot/b(i)
               endif
            enddo

         enddo
      enddo

      do i=1,n
         a(i,i)=1.0/b(i)
         do j=i+1,n
            tot=0.0
            do k=i,j-1
               tot=tot-a(j,k)*a(k,i)
            enddo
            a(j,i)=tot/b(j)
         enddo
      enddo

      do j=1,n
         do i=j,n
            tot=0.0
            do k=i,n
               tot=tot+a(k,i)*a(k,j)
            enddo
            if(i.eq.j) then
               a(i,i)=tot
            else
               a(i,j)=tot
               a(j,i)=a(i,j)
            endif

         enddo
      enddo

end subroutine PDmatinv
end module pwLKH
