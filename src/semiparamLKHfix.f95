!######################################################################

    module semiparamLKHfix
    use, intrinsic :: iso_c_binding
    implicit none
    public :: ORMLECOLfix, orMLErowfix

    contains

!######################################################################
!--|------------------------------------------|
!  | Direct Newton Raphson method for finding |
!  | the maximizer in the linear OR model     |
!  |   Reference point prespecified           |
!--|------------------------------------------|
!    gam: odds ratio estimator
!    estv: variance estimate
!---------------------------------------------------------

!-------------------------------------------
! Parameter gam is vectorized by column
!-------------------------------------------
Subroutine ORMLECOLfix(y,x,n,np,nq,fix,gam,estv,loglkh,niter,eps,converge) bind(C, name = "ormlecolfix_")

      implicit none
      !integer, parameter:: dp = selected_real_kind(15, 307)

      integer (C_INT) n,np,nq,nc,niter,converge! integer reference point
      integer (C_INT) fix(np,nq)
      real (C_DOUBLE) y(n,np),x(n,nq)
      real (C_DOUBLE) gam(np*nq),estv(np*nq,np*nq),loglkh
              ! row-wise vectorization

      integer (C_INT) index(n),cat(n+1) ! indicator where an odds ratio parameter is set to zero
      real (C_DOUBLE) eps,delta,maxcp
      integer (C_INT) i,iarep,j,jj,k,ki,kj,kk,kkjj,kl,l,m,na

      real (C_DOUBLE) TINY
      parameter(TINY=1e-6)

      !integer,allocatable:: Model(:) ! indicator of nonzero components in the parameter matrix
      integer,allocatable::sfix(:) ! straightened structure
      real (C_DOUBLE),allocatable:: dint(:),cp(:,:)!conditional probability
      real (C_DOUBLE),allocatable:: ce(:,:),cov(:,:,:)!,vtg(np*nq,np*nq) !conditional mean and variance
      real (C_DOUBLE),allocatable:: der(:),der2(:,:) !first and second derivatives
      real (C_DOUBLE),allocatable::  para(:),tp(:)  !total parameter vector

      call orderANDcount(y,n,np,index,cat(1:(nc+1)),nc)
      !write(*,*)"# of cat",nc
      na=np*nq
      m=np*nq+nc-1
      allocate(dint(n),cp(n,nc),para(nc))
      allocate(ce(n,np),cov(n,np,np),der(m),der2(m,m),tp(m),sfix(m))

      do k=1,nc
         para(k)=log(1.0*(cat(k+1)-cat(k))/(cat(nc+1)-cat(nc))) !nonparametric parameter
      enddo
      !do k=1,na !initial values are important to keep
      !   gam(k)=0.0  !odds ratio parameter
      !enddo
      tp=0

      do k=1,nc-1
        sfix(k)=1 !for baseline
      enddo
      do k=1,nq
        do kk=1,np
          ki=(k-1)*np+kk
          sfix(ki+nc-1)=fix(kk,k)
        enddo
      enddo

    Iteration:DO IArep=1,niter  ! overall iteration

      do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
         dint(i)=0.0
         do k=1,np
            ce(i,k)=0.0
            do l=1,np
               cov(i,k,l)=0.0
            enddo
         enddo

         do j=1,nc    ! compute denominator integral
            cp(i,j)=0
            do k=1,nq
               do l=1,np
                  !if(fix(l,k)/=0) then
                    kl=(k-1)*np+l
                    cp(i,j)=cp(i,j)+gam(kl)*y(index(cat(j)),l)*x(i,k)
                  !endif
               enddo
            enddo
         enddo

         maxcp=maxval(cp(i,:)+para) ! subtract maximum cp value for each subject to achieve numeric stability
         do j=1,nc
            cp(i,j)=exp(cp(i,j)+para(j)-maxcp)
            dint(i)=dint(i)+cp(i,j) !for denominator

            do k=1,np  ! for conditional probability, mean, and variance
               ce(i,k)=ce(i,k)+y(index(cat(j)),k)*cp(i,j)
               do l=1,np
                  cov(i,k,l)=cov(i,k,l)+y(index(cat(j)),k)*y(index(cat(j)),l)*cp(i,j)
               enddo
            enddo

         enddo

         do j=1,nc  ! CONDITIONAL PROBABILITY
            cp(i,j)=cp(i,j)/dint(i)
         enddo

         do k=1,np ! complete the conditional expectation calculation
            ce(i,k)=ce(i,k)/dint(i)
            do l=1,k
               cov(i,k,l)=cov(i,k,l)/dint(i)-ce(i,k)*ce(i,l)
               cov(i,l,k)=cov(i,k,l)
            enddo
         enddo

      enddo

     ! Find the first and second derivatives for nonparametric part

      do j=1,m
         if(j<nc) then
            der(j)=cat(j+1)-cat(j)
         else
            der(j)=0.0
         endif
         do k=1,j
            der2(j,k)=0.0
            der2(k,j)=0.0
         enddo
      enddo

      subject:do i=1,n

         do j=1,nc-1
            der(j)=der(j)-cp(i,j)
            do l=1,j-1
               der2(j,l)=der2(j,l)+cp(i,j)*cp(i,l)
            enddo
            der2(j,j)=der2(j,j)-cp(i,j)*(1-cp(i,j))
         enddo

      !
      ! Find the first and second deivatives for parametric part
      !

         do k=1,nq
            do j=1,np

              if(fix(j,k)/=0) then
                kj=nc-1+(k-1)*np+j
                der(kj)=der(kj)+(y(i,j)-ce(i,j))*x(i,k)

                do kk=1,k
                  do jj=1,np
                    if(fix(jj,kk)/=0)then
                      kkjj=nc-1+(kk-1)*np+jj
                      if(kkjj<=kj) then
                        der2(kj,kkjj)=der2(kj,kkjj)-cov(i,j,jj)*x(i,k)*x(i,kk)
                      endif
                    endif
                  enddo
                enddo

              endif

            enddo
         enddo

      !
      ! Find the joint second derivative of nonparametric
      !       and parametric parts
      !
         do k=1,nq
            do j=1,np
              if(fix(j,k)/=0) then
                kj=nc-1+(k-1)*np+j
                do l=1,nc-1
                  der2(kj,l)=der2(kj,l)-(y(index(cat(l)),j)-ce(i,j))*cp(i,l)*x(i,k)
                enddo
              endif
            enddo
         enddo

      enddo subject  ! sum over subjects

      ! reflecting the second derivatives
      do j=1,m   ! m=nc-1+na+nb; na=np*nq; nb=np*nr
        do l=1,j-1
            der2(j,l)=-der2(j,l)
            der2(l,j)=der2(j,l)
          enddo
          der2(j,j)=-der2(j,j)
      enddo

!!!! Find the inverse of the second derivative matrix

      Call PDmatinvfix(der2,sfix,m)

!!!! parameter updating

      delta=0
      do j=1,m

         tp(j)=0
         do l=1,m
           tp(j)=tp(j)+der2(j,l)*der(l)
         enddo
         delta=delta+abs(tp(j))

      enddo

      if(iarep<10) tp=tp/4

      do j=1,nc-1
         para(j)=para(j)+tp(j)
      enddo
      do k=1,nq
        do j=1,np
          if(fix(j,k)/=0) then
            kj=(k-1)*np+j
            gam(kj)=gam(kj)+tp(nc-1+kj)
          endif
        enddo
      enddo

!      write(*,*) IArep,delta
!      write(*,*) (gam(j),j=1,na)
!      write(*,*) (para(j),j=1,nc)
!      pause
!!!! check convergence

      if(delta<eps) Then

         ! convergent
         ! (a). Find the variance for the gamma estimate
        converge=1
        call Likelihoodcol(y,x,n,np,nq,gam,para,index,cat,nc,loglkh)

        ! subset the variance matrix
        do k=1,nq
          do j=1,np

            if(fix(j,k)/=0) then
              kj=(k-1)*np+j
              do kk=1,nq
                do jj=1,np
                  if(fix(jj,kk)/=0) then
                    kkjj=(kk-1)*np+jj
                    estv(kj,kkjj)=der2(kj+nc-1,kkjj+nc-1)
                  endif
                enddo
              enddo
            endif

          enddo
        enddo

        exit iteration

      endif

   enddo iteration


   Deallocate(dint,cp,para,ce,cov,der,der2,tp,sfix)

 end subroutine ORMLECOLfix

!---------------------------------------------------------
!   Compute log-likelihood at given parameter values(col)
!---------------------------------------------------------
subroutine Likelihoodcol(y,x,n,np,nq,gam,para,index,cat,nc,loglkh)

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc
   integer index(n),cat(nc+1) ! nc---total categories of y
                              ! index(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc)

   integer i,j,k,l
   real(kind=dp) dint,cp(nc),loglkh
   real(kind=dp) maxcp

   loglkh=0
   do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
      do j=1,nc    ! compute denominator integral
         cp(j)=0
         do k=1,nq
            do l=1,np
               cp(j)=cp(j)+gam((k-1)*np+l)*y(index(cat(j)),l)*x(i,k)
            enddo
         enddo
      enddo

      maxcp=maxval(cp+para) ! subtract maximum cp value for each subject to achieve numeric stability
      dint=0
      do j=1,nc
         dint=dint+exp(cp(j)+para(j)-maxcp)*(cat(j+1)-cat(j))/n !for denominator
      enddo
      do k=1,nq
         do l=1,np
            loglkh=loglkh+gam((k-1)*np+l)*y(i,l)*x(i,k)
         enddo
      enddo
      loglkh=loglkh-maxcp-log(dint)
   enddo
   do j=1,nc  ! add nonparametric part in the numerator of the likelihood
      loglkh=loglkh+(cat(j+1)-cat(j))*(para(j)+log(1.0*(cat(j+1)-cat(j))/n))
   enddo

end subroutine Likelihoodcol


!-------------------------------------------
! Parameter gam is vectorized by row.
!-------------------------------------------
Subroutine orMLErowfix(y,x,n,np,nq,fix,gam,estv,loglkh,niter,eps,converge) bind(C, name = "ormlerowfix_")


  implicit none
  !integer, parameter:: dp = selected_real_kind(15, 307)

      integer (C_INT) converge !=1 if converges, 0 otherwise.
      integer (C_INT) n,np,nq,nc,niter! integer reference point
      integer (C_INT)fix(np,nq)
      real (C_DOUBLE) y(n,np),x(n,nq)
      real (C_DOUBLE) gam(np*nq),estv(np*nq,np*nq),loglkh
              ! row-wise vectorization

      integer (C_INT) index(n),cat(n+1),Rindex(n) ! indicator where an odds ratio parameter is set to zero
      real (C_DOUBLE) eps,delta,maxcp
      integer (C_INT) i,iarep,j,jj,k,ki,kj,kk,kkjj,kl,l,m,na

      !integer,allocatable:: Model(:) ! indicator of nonzero components in the parameter matrix
      integer (C_INT),allocatable::sfix(:) ! straightened structure
      real (C_DOUBLE),allocatable:: dint(:),cp(:,:)!conditional probability
      real (C_DOUBLE),allocatable:: ce(:,:),cov(:,:,:)!,vtg(np*nq,np*nq) !conditional mean and variance
      real (C_DOUBLE),allocatable:: der(:),der2(:,:) !first and second derivatives
      real (C_DOUBLE),allocatable::  para(:),tp(:)  !total parameter vector

      call orderANDcount(y,n,np,index,cat,nc)

      converge=0
      na=np*nq
      m=np*nq+nc-1
      allocate(dint(n),cp(n,nc),para(nc))
      allocate(ce(n,np),cov(n,np,np),der(m),der2(m,m),tp(m),sfix(m))

      do k=1,nc
         para(k)=log(1.0*(cat(k+1)-cat(k))/(cat(nc+1)-cat(nc))) !nonparametric parameter
      enddo
      !do k=1,na ! initial values are important to keep
      !   gam(k)=0.0  !odds ratio parameter
      !enddo
      tp=0

      do k=1,nc-1
        sfix(k)=1 !for baseline
      enddo
      do k=1,np
        do kk=1,nq
          ki=(k-1)*nq+kk
          sfix(ki+nc-1)=fix(k,kk)
        enddo
      enddo

    Iteration:DO IArep=1,niter  ! overall iteration

      do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
         dint(i)=0.0
         do k=1,np
            ce(i,k)=0.0
            do l=1,np
               cov(i,k,l)=0.0
            enddo
         enddo

         do j=1,nc    ! compute denominator integral
            cp(i,j)=0
            do k=1,np
               do l=1,nq
                  !if(fix(k,l)/=0) then
                    kl=(k-1)*nq+l
                    cp(i,j)=cp(i,j)+gam(kl)*y(index(cat(j)),k)*x(i,l)
                  !endif
               enddo
            enddo
         enddo

         maxcp=maxval(cp(i,:)+para) ! subtract maximum cp value for each subject to achieve numeric stability
         do j=1,nc
            cp(i,j)=exp(cp(i,j)+para(j)-maxcp)
            dint(i)=dint(i)+cp(i,j) !for denominator

            do k=1,np  ! for conditional probability, mean, and variance
               ce(i,k)=ce(i,k)+y(index(cat(j)),k)*cp(i,j)
               do l=1,np
                  cov(i,k,l)=cov(i,k,l)+y(index(cat(j)),k)*y(index(cat(j)),l)*cp(i,j)
               enddo
            enddo

         enddo

         do j=1,nc  ! CONDITIONAL PROBABILITY
            cp(i,j)=cp(i,j)/dint(i)
         enddo

         do k=1,np ! complete the conditional expectation calculation
            ce(i,k)=ce(i,k)/dint(i)
            do l=1,k
               cov(i,k,l)=cov(i,k,l)/dint(i)-ce(i,k)*ce(i,l)
               cov(i,l,k)=cov(i,k,l)
            enddo
         enddo

      enddo

     ! Find the first and second derivatives for nonparametric part

      do j=1,m
         if(j.lt.nc) then
            der(j)=cat(j+1)-cat(j)
         else
            der(j)=0.0
         endif
         do k=1,j
            der2(j,k)=0.0
            der2(k,j)=0.0
         enddo
      enddo

      subject:do i=1,n

         do j=1,nc-1
            der(j)=der(j)-cp(i,j)
            do l=1,j-1
               der2(j,l)=der2(j,l)+cp(i,j)*cp(i,l)
            enddo
            der2(j,j)=der2(j,j)-cp(i,j)*(1-cp(i,j))
         enddo

      !
      ! Find the first and second deivatives for parametric part
      !

         do k=1,np
            do j=1,nq

              if(fix(k,j)/=0) then
                kj=nc-1+(k-1)*nq+j
                der(kj)=der(kj)+(y(i,k)-ce(i,k))*x(i,j)
                do kk=1,k
                  do jj=1,nq
                    if(fix(kk,jj)/=0) then
                      kkjj=nc-1+(kk-1)*nq+jj
                      if(kkjj<=kj) then
                        der2(kj,kkjj)=der2(kj,kkjj)-cov(i,k,kk)*x(i,j)*x(i,jj)
                      endif
                    endif
                  enddo
                enddo
              endif

            enddo
         enddo

      !
      ! Find the joint second derivative of nonparametric
      !       and parametric parts
      !
         do k=1,np
            do j=1,nq

              if(fix(k,j)/=0) then
                kj=nc-1+(k-1)*nq+j
                do l=1,nc-1
                  der2(kj,l)=der2(kj,l)-(y(index(cat(l)),k)-ce(i,k))*cp(i,l)*x(i,j)
                enddo
              endif

            enddo
         enddo

      enddo subject  ! sum over subjects

      ! reflecting the second derivatives
      do j=1,m   ! m=nc-1+na+nb; na=np*nq; nb=np*nr
        do l=1,j-1
            der2(j,l)=-der2(j,l)
            der2(l,j)=der2(j,l)
          enddo
          der2(j,j)=-der2(j,j)
      enddo

!!!! Find the inverse of the second derivative matrix

      Call PDmatinvfix(der2,sfix,m)

!!!! parameter updating

      delta=0
      do j=1,m

         tp(j)=0
         do l=1,m
           tp(j)=tp(j)+der2(j,l)*der(l)
         enddo
         delta=delta+abs(tp(j))

      enddo

   if(iarep<10) tp=tp/4

      do j=1,nc-1
         para(j)=para(j)+tp(j)
      enddo
      do k=1,np
        do j=1,nq
          if(fix(k,j)/=0) then
            kj=(k-1)*nq+j
            gam(kj)=gam(kj)+tp(nc-1+kj)
          endif
        enddo
      enddo

      if(delta.lt.eps) Then

        converge=1
         ! convergent
                 !compute loglikelihood
        call Likelihoodrow(y,x,n,np,nq,gam,para,index,cat,nc,loglkh)

         ! subset the variance matrix
        do k=1,np
          do j=1,nq

            if(fix(k,j)/=0) then
              kj=(k-1)*nq+j
              do kk=1,np
                do jj=1,nq
                  if(fix(kk,jj)/=0) then
                    kkjj=(kk-1)*nq+jj
                    estv(kj,kkjj)=der2(kj+nc-1,kkjj+nc-1)
                  endif
                enddo
              enddo
            endif

          enddo
        enddo

        exit iteration
      endif

   enddo iteration


   Deallocate(dint,cp,para,ce,cov,der,der2,tp,sfix)

 end subroutine orMLErowfix

 !---------------------------------------------------------
!   Compute log-likelihood at given parameter values(row)
!---------------------------------------------------------
subroutine Likelihoodrow(y,x,n,np,nq,gam,para,index,cat,nc,loglkh)

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc
   integer index(n),cat(nc+1) ! nc---total categories of y
                              ! index(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc)

   integer i,j,k,l
   real(kind=dp) dint,cp(nc),loglkh
   real(kind=dp) maxcp

   loglkh=0
   do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
      do j=1,nc    ! compute denominator integral
         cp(j)=0
         do k=1,np
            do l=1,nq
               cp(j)=cp(j)+gam((k-1)*nq+l)*y(index(cat(j)),k)*x(i,l)
            enddo
         enddo
      enddo

      maxcp=maxval(cp+para) ! subtract maximum cp value for each subject to achieve numeric stability
      dint=0
      do j=1,nc
         dint=dint+exp(cp(j)+para(j)-maxcp)*(cat(j+1)-cat(j))/n !for denominator
      enddo
      do k=1,np
         do l=1,nq
            loglkh=loglkh+gam((k-1)*nq+l)*y(i,k)*x(i,l)
         enddo
      enddo
      loglkh=loglkh-maxcp-log(dint)
   enddo
   do j=1,nc  ! add nonparametric part in the numerator of the likelihood
      loglkh=loglkh+(cat(j+1)-cat(j))*(para(j)+log(1.0*(cat(j+1)-cat(j))/n))
   enddo

end subroutine Likelihoodrow


!----------------------------------------------
!  Order data and get location or categories
!----------------------------------------------
 subroutine orderANDcount(x,n,p,index,cat,nc)

  implicit none
  integer, parameter:: dp = selected_real_kind(15, 307)

      real TINY
      parameter(TINY=10e-6)
      integer n,p,nc
      integer index(n),cat(n+1)
      real(kind=dp) x(n,p)

      integer itemp,i
      real rtemp

      call order(x,n,p,index)

      rtemp=x(index(1),p)  ! collect the numbers of categories
      itemp=1
      cat(1)=1
      do i=2,n
         if(abs(x(index(i),p)-rtemp).GT.TINY) then
            rtemp=x(index(i),p)
            itemp=itemp+1
            cat(itemp)=i   ! The beginning index.
         endif
      enddo
      nc=itemp
      cat(nc+1)=n+1
     ! write(*,*)(cat(i), i=1,nc+1)

   end subroutine orderANDcount

   !---|---------------------|
   !   | Sorting subroutine  |
   !---|---------------------|---------------------------
   !    x(n,p): input sequence
   !    index(n): output of address order of sequence
   !-----------------------------------------------------
   subroutine order(x,n,p,index)


  implicit none
  integer, parameter:: dp = selected_real_kind(15, 307)


      integer n,p,index(n)
      real(kind=dp) x(n,p)

      real TINY
      parameter(TINY=10e-6)
      integer low,upp,indtemp,i,j,ii,l
      real temp

      do i=1,n
         index(i)=i
      enddo

      do i=1,(n-1)
         temp=x(index(i),1)
         do j=(i+1),n
            if(temp.GT.x(index(j),1))then
               temp=x(index(j),1)
               indtemp=index(i)
               index(i)=index(j)
               index(j)=indtemp
            else
            endif
         enddo
      enddo


      do j=2,p

         low=1
         temp=x(index(1),j-1)
         upp=low
         do i=2,n
            if(abs(x(index(i),j-1)-temp).gt.TINY) then
               if(upp.gt.low) then
               do ii=low,upp
                  temp=x(index(ii),j)
                  do l=(ii+1),upp
                     if(temp.gt.x(index(l),j)) then
                        temp=x(index(l),j)
                        indtemp=index(ii)
                        index(ii)=index(l)
                        index(l)=indtemp
                     endif
                  enddo
               enddo
               endif
               temp=x(index(i),j-1)  !start next ordering
               low=upp+1
               upp=low
            elseif(i.eq.n) then
               if(n.gt.low) then
               do ii=low,n
                  temp=x(index(ii),j)
                  do l=(ii+1),n
                     if(temp.gt.x(index(l),j)) then
                        temp=x(index(l),j)
                        indtemp=index(ii)
                        index(ii)=index(l)
                        index(l)=indtemp
                     endif
                  enddo
               enddo
               endif
            else
               upp=upp+1 ! counting within a category
            endif
         enddo
      enddo

   end subroutine order


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
end module semiparamLKHfix
