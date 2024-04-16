!######################################################################

    module semiparamPENLKH
    use, intrinsic :: iso_c_binding
    implicit none
    public :: NetworkSelectbySP

    contains

!######################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Network selection by
!!    penalized semiparametric likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------
!            Network selection
!--------------------------------------------------------------------------
! Input
!   dat(n,node)==> input data
!   ng ==> number of groups the nodes are divided.
!   group(ng) ==> number of nodes in each group
!   nlam ==> number of different penalty values
!   lambda(nlma) ==> different penalty values
! Output
!   Selnet==> selection outcomes Sel(nlam,node,node)
!-------------------------------------------------------------------------
 subroutine NetworkSelectbySP(dat,n,node,group,ng,lam,nlam,niter,eps,selnet) bind(C, name = "networkselectbysp_")

   !integer, parameter:: dp = selected_real_kind(15, 307)

   integer (C_INT) n,node,ng   ! ng ==> total number of clusters
   integer (C_INT) group(ng) ! size of each cluster
   real (C_DOUBLE) dat(n,node) !input data
   integer (C_INT) cat(n+1),index(n),nc
   integer (C_INT) np,nq

   integer (C_INT) selnet(node,node,nlam)

   integer (C_INT) nlam,niter  !nlam is the number of lambda values to search
   real (C_DOUBLE) eps

   real (C_DOUBLE) lam(nlam)

   integer (C_INT) j,k,ncurrent
   integer (C_INT) init

   real (C_DOUBLE),allocatable:: datswp(:,:) ! used for data swap: datswp(n,np)
   real (C_DOUBLE),allocatable:: para(:),gam(:),mgam(:,:)

  !!!! NETWORK DATA ANALYSIS !!!!


   selnet=0
   ncurrent=0 ! current location of the cluster

   do k=1,ng

      np=group(k)
      nq=node-np

      ALLOCATE(datswp(n,np))

      datswp=dat(:,(ncurrent+1):(ncurrent+np))
      dat(:,(ncurrent+1):(ncurrent+np))=dat(:,1:np)
      dat(:,1:np)=datswp

      ! order and count the categories for outcome data
      call orderANDcount(dat(:,1:np),n,np,index,cat,nc)

      ALLOCATE(para(nc),gam(np*nq),mgam(np,nq))

      gam=0.0
      para=0.0
      init=1
      do j=1,nlam
         ! selection by penalized semiparametric likelihood
         call ORCordSel(dat(:,1:np),dat(:,(np+1):node),n,np,nq,gam,para,index,cat,nc, &
                        selnet((ncurrent+1):(ncurrent+np),1:(node-np),j),niter,eps,lam(j),init) ! 0 mean parameter has already been initialized
      enddo
     !restore original data.
     !NOTE the order of the two statements are important and should not be changed
      dat(:,1:np)=dat(:,(ncurrent+1):(ncurrent+np))
      dat(:,(ncurrent+1):(ncurrent+np))=datswp

      DEALLOCATE(datswp,para,gam,mgam)

      ncurrent=ncurrent+np ! move to next cluster

   enddo

  ! Adjust positions of selected model index due to data sweep distortion
    do j=1,nlam !adjust position of select model index for data sweep distortion
     ncurrent=0
     do k=1,ng ! adjust positions
       np=group(k)
       if(k<ng) then
   selnet((ncurrent+1):(ncurrent+np),(ncurrent+np+1):node,j)=selnet((ncurrent+1):(ncurrent+np),(ncurrent+1):(node-np),j)
       endif
       selnet((ncurrent+1):(ncurrent+np),(np+1):(ncurrent+np),j)=selnet((ncurrent+1):(ncurrent+np),1:ncurrent,j)
       selnet((ncurrent+1):(ncurrent+np),1:np,j)=selnet((ncurrent+1):(ncurrent+np),(ncurrent+1):(ncurrent+np),j)
       selnet((ncurrent+1):(ncurrent+np),(ncurrent+1):(ncurrent+np),j)=0
       ncurrent=ncurrent+np
     enddo
     selnet(:,:,j)=selnet(:,:,j)+transpose(selnet(:,:,j))
   enddo

end subroutine NetworkSelectbySP

!--|------------------------------------------|
!  | Variable selection in SPORM              |
!--|------------------------------------------|
!    gam: odds ratio estimator
!    para: baseline parameter
!       sel --- the selection indicator
!       lambda --- tuning parameter
!-----------------------------------------------

Subroutine ORCordSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,lambda,init)

    integer, parameter:: dp = selected_real_kind(15, 307)

      integer n,np,nq,nc,niter,init ! integer reference point
      integer indexs(n),cat(n+1),sel(np,nq)

      real(kind=dp) y(n,np),x(n,nq)
      real(kind=dp) eps,lambda
      real(kind=dp) gam(np*nq),para(nc)

      integer nstep,nchange ! nchange is convergence monitor

      integer irep,k

      if(init==1)then
        do k=1,nc !nonparametric parameter
           para(k)=log(1.0*(cat(k+1)-cat(k))/n)
        enddo
        gam=0.0   ! OR parameters
        sel=0   ! OR selection indicator
      else
        !assume gam is given, estimate baseline first
        !call BaseEstbyApprxSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)
        call BaseEstbyNRSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)
      endif

      iteration:DO irep=1,niter  ! overall iteration

         !write(*,*)irep,irep
         !pause
         !1. select variable

         call ORselect(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nchange,lambda)
         !call MORselect(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nchange,lambda)

         !2. update baseline
         !if(irep<(200/lambda)) then
         !if(irep<10) then
           call BaseEstbyApprxSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)
         !else
         !  call BaseEstbyNRSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)
         !endif

!write(*,*)irep,nchange,nstep,sum(sel)

!         if(irep>(200/lambda).and.nchange==0.and.nstep==1)then
         if(irep>10.and.nchange==0.and.nstep==1)then
           ! write(*,*)irep,nchange,nstep,sum(sel)
            exit iteration
         !elseif(irep==niter) then
         !   write(*,*)irep,nchange,nstep,sum(sel),"unconvergent"
         endif

      enddo iteration

  !    do i=1,np
  !       do k=1,nq
  !          ik=(i-1)*nq+k
  !          if(sel(i,k)==1) then
  !             write(*,11) ik, gam(ik)
  !             11 format(I10,1x,20(es10.2,2x))
  !          endif
  !       enddo
  !    enddo

end subroutine ORCordSel


!-------------------------------------------------------
!    Odds Ratio model selection by l1-penalization
!-------------------------------------------------------
Subroutine ORselect(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nchange,lambda)

   !Use compORLoglikelihood,only:loglkh

   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc,niter,nchange ! niter-- maximum nuber of iteration; nstep--actual number of steps used
   integer indexs(n),cat(nc+1) ! nc---total categories of y
                              ! indexs(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   real(kind=dp) eps,lambda         ! lambda -- penalty
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc)
   integer sel(np,nq)      ! indicator for odds ratio parameter selection

   integer i,j,k,l,ij,iv,jv,irep
   integer selold ! old selection indicator

   real(kind=dp) delta
   real(kind=dp) maxcp,dint,cp(nc),ce,cov,der,der2

   ! one round of selection

   nchange=0 ! count selected variable change

   Paramy:do iv=1,np ! loop through y
      Paramx:do jv=1,nq   ! loop through x
        !compute the first derivative first,
        !    1. if the absolute value is less that lambda, then set the corresponding theta zero and go to next parameter
        !    2. if the absolute value is greater than lambda, then Use one-step Harley approximation x=x0-[f(x0)/f'(x0)]/{1-0.5*f"(x0)*f(x0)/[f'(x0)]^2}.

        ij=(iv-1)*nq+jv  ! current parameter location in gam

        selold=sel(iv,jv) ! keep old selection indicator
        delta=10 ! set initial delta big
        iteration: do irep=1,niter
          der=0 ! compute the first derivative
          der2=0
          do i=1,n ! initialize conditional probability, expectation, and covariance of Y.

            do j=1,nc    ! compute denominator integral
              cp(j)=0
              do k=1,np
                do l=1,nq
                  if(sel(k,l)==1)then
                    cp(j)=cp(j)+gam((k-1)*nq+l)*y(indexs(cat(j)),k)*x(i,l)
                  endif
                enddo
              enddo
            enddo

            dint=0.0
            ce=0.0
            cov=0.0
            maxcp=maxval(cp+para)
            do j=1,nc
              cp(j)=exp(cp(j)+para(j)-maxcp)

              dint=dint+cp(j) !for denominator

             ! for conditional mean and variance
              ce=ce+y(indexs(cat(j)),iv)*cp(j)
              cov=cov+y(indexs(cat(j)),iv)*y(indexs(cat(j)),iv)*cp(j)
            enddo

            ! complete the conditional expectation calculation
            ce=ce/dint
            cov=cov/dint-ce*ce

             ! Find the first and second deivatives for parametric part.
            der=der-(y(i,iv)-ce)*x(i,jv)
            der2=der2+cov*x(i,jv)*x(i,jv)
          enddo

          if(gam(ij)>(der+lambda)/der2) then
            delta=(der+lambda)/der2
            gam(ij)=gam(ij)-delta
            sel(iv,jv)=1
          elseif(gam(ij)<(der-lambda)/der2) then
            delta=(der-lambda)/der2
            gam(ij)=gam(ij)-delta
            sel(iv,jv)=1
          else
            delta=abs(gam(ij))
            gam(ij)=0
            sel(iv,jv)=0

          endif

      !    if(sel(iv,jv)==1) then
      !       write(*,*) ij,delta,gam(ij)
      !    endif

          if(abs(delta)<eps)then
            exit iteration
          endif

        enddo iteration

        nchange=nchange+abs(selold-sel(iv,jv)) ! count selection changes

         !write(*,*) 'Number of variables selected: ', sum(sel)

     enddo Paramx
   enddo Paramy

end subroutine ORselect

!-------------------------------------------------------
!    Odds Ratio model selection by l1-penalization
!    using multivariate quadratic programming
!-------------------------------------------------------
Subroutine MORselect(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nchange,lambda)

   implicit none

   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc,niter,nchange ! niter-- maximum nuber of iteration; nstep--actual number of steps used
   integer indexs(n),cat(nc+1) ! nc---total categories of y
                              ! indexs(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   real(kind=dp) lambda,eps         ! lambda -- penalty
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc)
   integer sel(np,nq),vsel(np*nq)      ! indicator for odds ratio parameter selection


   integer i,j,k,l,j1,j2,jj,k1,k2,kk

   real(kind=dp) maxcp,dint,cp(nc),ce(np),cov(np,np),der(np*nq),der2(np*nq,np*nq)

   der=0 ! compute the first derivative
   der2=0
   do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
     dint=0.0
     ce=0.0
     cov=0.0

     do j=1,nc    ! compute denominator integral
       cp(j)=0
       do k=1,np
         do l=1,nq
           if(sel(k,l)==1)then
             cp(j)=cp(j)+gam((k-1)*nq+l)*y(indexs(cat(j)),k)*x(i,l)
           endif
         enddo
       enddo
     enddo

     maxcp=maxval(cp+para)
     do j=1,nc
       cp(j)=exp(cp(j)+para(j)-maxcp)
       dint=dint+cp(j) !for denominator

      ! for conditional mean and variance
       do jj=1,np
         ce(jj)=ce(jj)+y(indexs(cat(j)),jj)*cp(j)
         do kk=1,np
           cov(jj,kk)=cov(jj,kk)+y(indexs(cat(j)),jj)*y(indexs(cat(j)),kk)*cp(j)
         enddo
       enddo

     enddo

     ! complete the conditional expectation calculation
     do jj=1,np
       ce(jj)=ce(jj)/dint
       do kk=1,np
         cov(jj,kk)=cov(jj,kk)/dint-ce(jj)*ce(kk)
       enddo
     enddo

     ! Find the first and second deivatives for parametric part.
     do j1=1,np
       do j2=1,nq
         jj=(j1-1)*nq+j2
         der(jj)=der(jj)-(y(i,j1)-ce(j1))*x(i,j2)

         do k1=1,np
           do k2=1,nq
             kk=(k1-1)*nq+k2
             der2(jj,kk)=der2(jj,kk)+cov(j1,k1)*x(i,j2)*x(i,k2)
           enddo
         enddo

       enddo
     enddo

   enddo  ! complete calculating the first and second derivatives

!write(*,*)'der=',der
!write(*,*)'der2=',der2
   ! Do parameter selection by quadratic programming
   call QuadP(der,der2,gam,vsel,np*nq,niter,lambda)
   nchange=sum(vsel)-sum(sel) ! number of variables selected
   if(nchange<0.and.sum(vsel)==0) then  !reset initial values
     para=0.0
     gam=0.0
   endif
   !convert vectorized selection matrix back to matrix

   do j=1,np
     do k=1,nq
       jj=(j-1)*nq+k
       sel(j,k)=vsel(jj)
     enddo
   enddo

!write(*,*)sum(sel),sum(vsel)

end subroutine MORselect

!------------------------------------------------------------------------
!  Quadratic programming for estimating theta parameters
!  Input: mu (mean vector) and sig (the variance matrix),
!         lambda (penalty), and theta0.
!  Output: Parameter estimator for theta
!  The problem: min_{theta} mu*(theta-theta0)
!                          +0.5*(theta-theta0)Sig(theta-theta0)
!  Given the current theta and current selected variables in SEL,
!       update both beta and SEL
!-------------------------------------------------------------------------
subroutine QuadP(mu,sig,theta,sel,q,niter,lambda)

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer q,sel(q),niter
   real(kind=dp) mu(q),sig(q,q),theta(q),lambda

   integer k,irep,jv,nselect,selold
   real(kind=dp) delta,theta0(q)

   theta0=theta ! record the initial values
   nselect=0  ! set newly selected to 0.

   do irep=1,niter

     do jv=1,q ! go through all the variables
       selold=sel(jv)

       delta=sig(jv,jv)*theta0(jv)-mu(jv)
       do k=1,q
         delta=delta-sig(jv,k)*(theta(k)-theta0(k))
       enddo
       delta=delta+sig(jv,jv)*(theta(jv)-theta0(jv)) ! remove the extra item should not be included

       if(delta>lambda) then
         theta(jv)=(delta-lambda)/sig(jv,jv)
         Sel(jv)=1
       elseif(delta<-lambda) then
         theta(jv)=(delta+lambda)/sig(jv,jv)
         Sel(jv)=1
       else
         theta(jv)=0
         Sel(jv)=0
       endif

       nselect=nselect+abs(sel(jv)-selold)
     enddo

     if(nselect==0) then
       exit
     endif

   enddo

end subroutine QuadP
!-------------------------------------------------------------------------------
!   Nonparametric estimation of baseline function by successive approximation
!      for fixed parameter values
!-------------------------------------------------------------------------------
subroutine BaseEstbyApprxSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)

    integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc,niter,nstep ! niter-- maximum nuber of iteration; nstep--actual number of steps used
   integer indexs(n),cat(nc+1) ! nc---total categories of y
                              ! indexs(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   integer sel(np,nq) ! selected odd ratio parameters
   real(kind=dp) eps,delta,maxcp
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc),paranew(nc)

   integer irep,i,j,k,l
   real(kind=dp) dint(n),cp(n,nc)


   sucApprx: do irep=1,niter

      do i=1,n ! initialize conditional probability, expectation, and covariance of Y.
         dint(i)=0.0
         do j=1,nc    ! compute denominator integral
            cp(i,j)=0
            do k=1,np
               do l=1,nq
                  if(sel(k,l)==1) then
                     cp(i,j)=cp(i,j)+gam((k-1)*nq+l)*y(indexs(cat(j)),k)*x(i,l)
                  endif
               enddo
            enddo
         enddo

         maxcp=maxval(cp(i,:)) ! subtract maximum cp value for each subject to achieve numeric stability
         do j=1,nc
            cp(i,j)=exp(cp(i,j)-maxcp)
            dint(i)=dint(i)+cp(i,j)*exp(para(j))*(cat(j+1)-cat(j))/n !for denominator
         enddo
      enddo

      do j=1,nc ! update nonparametric parameters
         paranew(j)=0
         do i=1,n
            paranew(j)=paranew(j)+cp(i,j)/dint(i)
         enddo
         paranew(j)=1.0/paranew(j)
      enddo

      !paranew=paranew/sum(paranew) !!! There is an error here
      paranew=log(1.0*paranew/sum(paranew))

      delta=sum(abs(paranew-para))

      !write(*,*)irep,delta

      if(delta<eps) then
        ! write(*,*) para
        nstep=irep
        exit sucApprx
      else
        para=paranew
      endif

   enddo sucApprx


end subroutine BaseEstbyApprxSel

!-------------------------------------------------------------------------------
!   Nonparametric estimation of baseline function by Newton-Raphson method
!      for fixed parameter values
!-------------------------------------------------------------------------------
subroutine BaseEstbyNRSel(y,x,n,np,nq,gam,para,indexs,cat,nc,sel,niter,eps,nstep)

    integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,np,nq,nc,niter,nstep! niter-- maximum nuber of iteration; nstep--actual number of steps used
   integer indexs(n),cat(nc+1) ! nc---total categories of y
                              ! indexs(i)---category of subject i
                              ! cat(j) starting location of category j in ordered subject list
   integer sel(np,nq) ! odds ratio parameters selected and used in the calculation
   real(kind=dp) eps
   real(kind=dp) y(n,np),x(n,nq)
   real(kind=dp) gam(np*nq),para(nc)

   integer irep,i,j,k,l
   real(kind=dp) delta
   real(kind=dp) dint(n),cp(n,nc)
   real(kind=dp) der(nc-1),der2(nc-1,nc-1),tp(nc-1)

   NR: DO irep=1,niter  ! overall iteration

      do i=1,n ! initialize conditional probability of Y.
         dint(i)=0.0
         do j=1,nc    ! compute denominator integral
            cp(i,j)=0
            do k=1,np
               do l=1,nq
                  if(sel(k,l)==1) then
                     cp(i,j)=cp(i,j)+gam((k-1)*nq+l)*y(indexs(cat(j)),k)*x(i,l)
                  endif
               enddo
            enddo
            cp(i,j)=exp(cp(i,j)+para(j))*(cat(j+1)-cat(j))/n
            dint(i)=dint(i)+cp(i,j)   !for denominator
         enddo

         do j=1,nc  ! CONDITIONAL PROBABILITY
            cp(i,j)=cp(i,j)/dint(i)
         enddo

      enddo

     !
     ! Find the first and second derivatives for nonparametric part
     !

      do j=1,nc-1
         der(j)=cat(j+1)-cat(j)
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

      enddo subject  ! sum over subjects

      ! reflecting the second derivatives
      do j=1,nc-1
         do l=1,j-1
            der2(j,l)=-der2(j,l)
            der2(l,j)=der2(j,l)
         enddo
         der2(j,j)=-der2(j,j)
      enddo

      !write(*,*)(der2(j,j),j=1,nc-1)
      !pause
!!!! Find the inverse of the second derivative matrix

      call chol(der2,tp,nc-1,nc-1)

!!!! parameter updating

      delta=0
      do j=1,nc-1
         tp(j)=0
         do l=1,nc-1
            tp(j)=tp(j)+der2(j,l)*der(l)
         enddo
         delta=delta+abs(tp(j))
      enddo

      if(irep<20) tp=tp/(20.0-irep)   !reduce stepsize initially to stablize iteration

      do j=1,nc-1
         para(j)=para(j)+tp(j)
      enddo


!     write(*,*) irep,delta
!      write(*,*) (para(j),j=1,nc)
!!!! check convergence

      if(delta<eps) Then
         nstep=irep
         exit NR
      endif

   enddo NR

end subroutine BaseEstbyNRSel



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! II. Data preprocessing and other relevant functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!1. sorting data
   !----------------------------------------------
   !  Order data and get location or categories
   !----------------------------------------------
   subroutine orderANDcount(x,n,p,index,cat,nc)

    integer, parameter:: dp = selected_real_kind(15, 307)


      integer n,p,nc
      integer index(n),cat(n+1)

      real(kind=dp) x(n,p)

      integer itemp,i
      real(kind=dp) rtemp

      call order(x,n,p,index)

      rtemp=x(index(1),p)  ! collect the numbers of categories
      itemp=1
      cat(1)=1
      do i=2,n
         if(abs(x(index(i),p)-rtemp).GT.1e-14) then
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

    integer, parameter:: dp = selected_real_kind(15, 307)


      integer n,p,index(n)
      real(kind=dp) x(n,p)

      integer low,upp,indtemp,i,j,ii,l
      real(kind=dp) temp

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
            if(abs(x(index(i),j-1)-temp).gt.1e-14) then
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

!!!2. matrix inversion

!--|----------------------------------------|
 !  | inversion of positive definite matrix  |
 !  | by cholesky decomposition              |
 !--|----------------------------------------|------------------
 !  | input a is the original matrix, output a is the inversion
 !  | the originalmatrix is not kept.
 !  | a(nxn) only a submatrix of mxm is inverted.
 !--------------------------------------------------------------
 subroutine chol(a,b,n,m)

    integer, parameter:: dp = selected_real_kind(15, 307)

      integer n,m
      real(kind=dp) a(n,n),b(n)

      integer i,j,k
      real(kind=dp) sum

      b(1)=sqrt(a(1,1))
      do j=2,m
        a(j,1)=a(j,1)/b(1)
      enddo

      do i=2,m
         do j=i,m
            sum=a(i,j)
            do k=i-1,1,-1
               sum=sum-a(i,k)*a(j,k)
               if(i.eq.j) then
                 ! if(sum.le.0) then
                 !    write(*,*)j,j
                 ! endif
                  b(i)=sqrt(sum)
               else
                  a(j,i)=sum/b(i)
               endif
            enddo

         enddo
      enddo

      do i=1,m
         a(i,i)=1.0/b(i)
         do j=i+1,m
            sum=0.0
            do k=i,j-1
               sum=sum-a(j,k)*a(k,i)
            enddo
            a(j,i)=sum/b(j)
         enddo
      enddo

      do j=1,m
         do i=j,m
            sum=0.0
            do k=i,m
               sum=sum+a(k,i)*a(k,j)
            enddo
            if(i.eq.j) then
               a(i,i)=sum
            else
               a(i,j)=sum
               a(j,i)=a(i,j)
            endif

         enddo
      enddo

 end subroutine chol
 end module semiparamPENLKH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of network selection by penalized
!! semiparametric likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
