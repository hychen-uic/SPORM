!######################################################################

    module pwLKHfix
    use, intrinsic :: iso_c_binding
    implicit none
    public :: pwMLECOLfix, pwMLErowfix

    contains

!######################################################################
!----------------------------------------------------------------------
!  pairwise likelihood approach
!  parameter gam is vectorized by column.
!    Fix provides the structure of the association to be computed.
!    NOTE: when the fix struct does not include parameter to estimate,
!          the method compute loglike value in this case.
!    Theta enables fixed parameter values (other than 0) for those 
!       not to be computed. 
!    If a theta component is not to be included in the model, 
!        it should be set to 0.
!----------------------------------------------------------------------

subroutine pwMLECOLfix(y,x,n,p,q,fix,theta,estv,loglkh,niter,eps,converge) bind(C, name = "pwmlecolfix_")

  implicit none
  !integer, parameter:: dp = selected_real_kind(15, 307)

  integer (C_INT) n,p,q,niter,converge
  integer (C_INT) fix(p,q)!provide information on which theta components are fixed at 0.
  integer (C_INT) sfix(p*q) !vectorized structure of fix(p,q)
  real (C_DOUBLE)  y(n,p),x(n,q)
  real (C_DOUBLE)  theta(p*q),estv(p*q,p*q),loglkh ! columnwise vectorization

  integer (C_INT) i,j,k,kk,l,ll,ki,li,irep
  real (C_DOUBLE)  temp(q,p),probtemp
  real (C_DOUBLE)  der(p*q),der2(p*q,p*q),delta(p*q),eps

  real (C_DOUBLE)  tder(n,p*q)

  ! compute sfix, the vectorized fix
  do k=1,q
    do kk=1,p
      ki=(k-1)*p+kk
      sfix(ki)=fix(kk,k)
    enddo
  enddo

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
            !if(fix(kk,k)/=0) then
              temp(k,kk)=(y(i,kk)-y(j,kk))*(x(j,k)-x(i,k))
              probtemp=probtemp+theta((k-1)*p+kk)*temp(k,kk)
            !endif
          enddo
        enddo
        if(probtemp<=0) then
          probtemp=exp(probtemp)/(1+exp(probtemp))
        else
          probtemp=1.0/(1+exp(-probtemp))
        endif

        do k=1,q
          do kk=1,p

            if(fix(kk,k)/=0) then
              ki=(k-1)*p+kk
              der(ki)=der(ki)+temp(k,kk)*probtemp

              do l=1,k
                do ll=1,p
                  if(fix(ll,l)/=0) then
                    li=(l-1)*p+ll
                    if(li<=ki) then
                      der2(ki,li)=der2(ki,li)+temp(k,kk)*temp(l,ll)*probtemp*(1.0-probtemp)
                    endif
                  endif
                enddo
              enddo
             endif

          enddo
        enddo

      enddo
    enddo

    do k=1,p*q
      do kk=1,k
        der2(kk,k)=der2(k,kk)
      enddo
    enddo

    if(sum(sfix)/=0) then  
      Call PDmatinvfix(der2,sfix,p*q)
      delta=matmul(der2,der)
    else ! The case with no parameter to estimate
      delta=0
    endif

    if(sum(abs(delta))<eps) then
      converge=1
      ! compute variance of the estimator using sandwich estimator
      ! and the pseudo-likelihood
      loglkh=0

      do i=1,n

        tder(i,:)=0
        do j=1,n
          probtemp=0
          do k=1,q
            do kk=1,p
              !if(fix(kk,k)/=0) then
                temp(k,kk)=(y(i,kk)-y(j,kk))*(x(i,k)-x(j,k))
                probtemp=probtemp+theta((k-1)*p+kk)*temp(k,kk)
              !endif
            enddo
          enddo
          if(probtemp<=0) then
            probtemp=exp(probtemp)/(1+exp(probtemp))
          else
            probtemp=1.0/(1+exp(-probtemp))
          endif
          loglkh=loglkh+log(probtemp)

          do k=1,q
            do kk=1,p
              if(fix(kk,k)/=0) then
                ki=(k-1)*p+kk
                tder(i,ki)=tder(i,ki)+temp(k,kk)*(1.0-probtemp)
              endif
            enddo
          enddo

        enddo
        tder=tder/(n-1)
      enddo

      do k=1,p*q
        if(sfix(k)/=0) then
          do kk=1,k
            if(sfix(kk)/=0)then
              estv(k,kk)=0
              do i=1,n
                estv(k,kk)=estv(k,kk)+tder(i,k)*tder(i,kk)
              enddo
              estv(k,kk)=estv(k,kk)*4/(n-1.0)
              estv(kk,k)=estv(k,kk)
            endif
          enddo
        endif
      enddo

      if(sum(sfix)/=0) then  ! exclude the case with no parameter to estimate
        der2=der2*n*(n-1)/2.0
        estv=matmul(der2,matmul(estv,der2)) ! sandwich estimate of variance
      endif 
        
      exit
    else
      theta=theta-delta
    endif

  enddo

end subroutine pwMLECOLfix
!----------------------------------------------------------------------
!  pairwise likelihood approach
!  parameter gam is vectorized by row.
!    Theta enables fixed parameter values (other than 0) for those 
!       not to be computed. 
!    If a theta component is not to be included in the model, 
!        it should be set to 0.
!-----------------------------------------------------------------------
subroutine pwMLErowfix(y,x,n,p,q,fix,theta,estv,loglkh,niter,eps,converge) bind(C, name = "pwmlerowfix_")

  implicit none
  !integer, parameter:: dp = selected_real_kind(15, 307)

  integer (C_INT) converge !=1 if converges, =0 otherwise.
  integer (C_INT) n,p,q,niter
  real (C_DOUBLE)  y(n,p),x(n,q)
  real (C_DOUBLE)  theta(p*q),estv(p*q,p*q),loglkh ! columnwise vectorization
  integer (C_INT) fix(p,q)!provide information on which theta components are fixed at 0.
  integer (C_INT) sfix(p*q) !vectorized structure of fix(p,q)

  integer (C_INT) i,j,k,kk,l,ll,ki,li,irep
  real (C_DOUBLE)  temp(p,q),probtemp
  real (C_DOUBLE)  der(p*q),der2(p*q,p*q),delta(p*q),eps

  real (C_DOUBLE)  tder(n,p*q)

! compute sfix, the vectorized fix
  do k=1,p
    do kk=1,q
      ki=(k-1)*q+kk
      sfix(ki)=fix(k,kk)
    enddo
  enddo

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
            !if(fix(k,kk)/=0) then
              temp(k,kk)=(y(i,k)-y(j,k))*(x(j,kk)-x(i,kk))
              probtemp=probtemp+theta((k-1)*q+kk)*temp(k,kk)
            !endif
          enddo
        enddo
        if(probtemp<=0) then
          probtemp=exp(probtemp)/(1+exp(probtemp))
        else
          probtemp=1.0/(1+exp(-probtemp))
        endif

        do k=1,p
          do kk=1,q
            if(fix(k,kk)/=0) then
              ki=(k-1)*q+kk
              der(ki)=der(ki)+temp(k,kk)*probtemp
              do l=1,k
                do ll=1,q
                  if(fix(l,ll)/=0) then
                    li=(l-1)*q+ll
                    if(li<=ki) then
                      der2(ki,li)=der2(ki,li)+temp(k,kk)*temp(l,ll)*probtemp*(1-probtemp)
                    endif
                  endif
                enddo
              enddo
            endif

          enddo
        enddo

      enddo
    enddo

    do k=1,p*q
      do kk=1,k
        der2(kk,k)=der2(k,kk)
      enddo
    enddo

    if(sum(sfix)/=0) then  
      Call PDmatinvfix(der2,sfix,p*q)
      delta=matmul(der2,der)
    else ! The case with no parameter to estimate
      delta=0
    endif
     
    if(sum(abs(delta))<eps) then
      converge=1
      ! compute variance of the estimator using sandwich estimator
      loglkh=0
      do i=1,n

        tder(i,:)=0
        do j=1,n
          probtemp=0
          do k=1,p
            do kk=1,q
              !if(fix(k,kk)/=0)then
                temp(k,kk)=(y(i,kk)-y(j,kk))*(x(i,k)-x(j,k))
                probtemp=probtemp+theta((k-1)*q+kk)*temp(k,kk)
              !endif
            enddo
          enddo
          if(probtemp<=0) then
            probtemp=exp(probtemp)/(1+exp(probtemp))
          else
            probtemp=1.0/(1+exp(-probtemp))
          endif
          loglkh=loglkh+log(probtemp)

          do k=1,p
            do kk=1,q
              if(fix(k,kk)/=0) then
                ki=(k-1)*q+kk
                tder(i,ki)=tder(i,ki)+temp(k,kk)*(1.0-probtemp)
              endif
            enddo
          enddo

        enddo
        tder=tder/(n-1)
      enddo

      do k=1,p*q
        do kk=1,k
          if(fix(k,kk)/=0) then
            estv(k,kk)=0
            do i=1,n
              estv(k,kk)=estv(k,kk)+tder(i,k)*tder(i,kk)
            enddo
            estv(k,kk)=estv(k,kk)*4/(n-1.0)

            estv(kk,k)=estv(k,kk)
          endif

        enddo
      enddo

      if(sum(sfix)/=0) then  ! exclude the case with no parameter to estimate
        der2=der2*n*(n-1)/2.0
        estv=matmul(der2,matmul(estv,der2)) ! sandwich estimate of variance
      endif 
 
      exit
    else
      theta=theta-delta
    endif

  enddo

end subroutine pwMLErowfix

!-------------------------------------------------------
!  Inverse a submatrix of a positive definite matrix
!    ind denotes the submatrix
!-------------------------------------------------------
subroutine PDmatinvfix(a,ind,n)

    implicit none
    integer, parameter:: dp = selected_real_kind(15, 307)

    integer n,ind(n)
    real(kind=dp) a(n,n)

    integer subn,subm,i,j
    integer,allocatable:: subind(:)
    real(kind=dp),allocatable:: suba(:,:)

    subn=0
    do i=1,n
      if(ind(i)/=0) then
        subn=subn+1
      endif
    enddo

    allocate(subind(subn),suba(subn,subn))

    subn=0
    do i=1,n
      if(ind(i)/=0) then
        subn=subn+1
        subind(subn)=i
        subm=subn-1 !discount the current index
        do j=i,n
          if(ind(j)/=0) then
            subm=subm+1
            suba(subn,subm)=a(i,j)
            suba(subm,subn)=a(i,j)
          endif
        enddo
      endif
    enddo

    call PDmatinv(suba,subn)

    do i=1,subn
      do j=i,subn
        a(subind(i),subind(j))=suba(i,j)
        a(subind(j),subind(i))=suba(i,j)
      enddo
    enddo

    Deallocate(subind,suba)

 end subroutine PDmatinvfix
!--------------------------------------------
! Inverse by Cholestry decomposition of
!       a positive definite matrix
!---------------------------------------------
subroutine PDmatinv(a,n)

      implicit none
      integer, parameter:: dp = selected_real_kind(15, 307)

      integer n,ind(n)

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
end module pwLKHfix
