!######################################################################

    module pwPENLKH
    use, intrinsic :: iso_c_binding
    implicit none
    public :: NetworkSelectbyPW,penpwlkh

    contains

!######################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Network analysis by pairwise likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------
!            Network selection by pairwise likelihood
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
 subroutine NetworkSelectbyPW(dat,n,node,group,ng,lam,nlam,niter,eps,selnet) bind(C, name = "networkselectbypw_")

   implicit none

   !integer, parameter:: dp = selected_real_kind(15, 307)

   integer (C_INT) n,node,nlam,niter,ng !nlam is the number of lambda values to search
   integer (C_INT) group(ng)
   real (C_DOUBLE) dat(n,node)
   real (C_DOUBLE) lam(nlam)

   real (C_DOUBLE) eps

   integer (C_INT) selnet(node,node,nlam)

   real (C_DOUBLE), allocatable:: gam(:),datswp(:,:)

   integer (C_INT) j,k,np,nq,ncurrent

  ! 2. Gaussian network model fit
   selnet=0
   ncurrent=0
   do k=1,ng  ! go through all nodes

      np=group(k)
      nq=node-np

      ALLOCATE(datswp(n,np),gam(np*nq))

      ! Switch two group of nodes as outcomes and covariates
      datswp=dat(:,(ncurrent+1):(ncurrent+np))
      dat(:,(ncurrent+1):(ncurrent+np))=dat(:,1:np)
      dat(:,1:np)=datswp

      !write(*,*) 'node=',k

      do j=1,nlam
      ! selection by penalized likelihood
        call penpwlkh(dat(:,1:np),dat(:,(np+1):node),n,np,nq,gam, &
                          selnet((ncurrent+1):(ncurrent+np),1:nq,j),niter,eps,lam(j))
      enddo

      !write(*,*) 'node=',k,' EBICnorm=',EBIC(:)

     !exchange location back
      dat(:,1:np)=dat(:,(ncurrent+1):(ncurrent+np))
      dat(:,(ncurrent+1):(ncurrent+np))=datswp

      DEallocate(datswp,gam)

      ncurrent=ncurrent+np

   enddo

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

 end subroutine NetworkSelectbyPW

!-----------------------------------------------------------------
!          Selection for multiple outcome model
!-----------------------------------------------------------------
! Penalized pairwise pseudo-likelihood for variable selection
!   with multiple outcomes (MO).
! Using l^1 penalty by coordinate descending approach
!-----------------------------------------------------------------
subroutine penpwlkh(y,x,n,p,q,theta,Sel,niter,eps,lam)

  implicit none

  integer, parameter:: dp = selected_real_kind(15, 307)

  integer n,p,q,niter
  integer Sel(p,q),selold
  real(kind=dp) lam,eps
  real(kind=dp) y(n,p),x(n,q)
  real(kind=dp) theta(p,q) ! columnwise vectorization

  integer i,j,k,kk,kv,lv,irep
  real(kind=dp) probtemp,thetaold
  real(kind=dp) der,der2,delta

  do kv=1,p
    do lv=1,q

      ITERATION:do irep=1,niter
        !write(*,*)irep,niter
        Selold=Sel(kv,lv)
        thetaold=theta(kv,lv)

        der=0
        der2=0
        do i=1,n
          do j=i+1,n

            probtemp=0
            do k=1,p
              do kk=1,q
                if(Sel(k,kk)==1)then
                  probtemp=probtemp+theta(k,kk)*(y(i,k)-y(j,k))*(x(i,kk)-x(j,kk))
                endif
              enddo
            enddo
            if(probtemp<=0) then
              probtemp=exp(probtemp)/(1+exp(probtemp))
            else
              probtemp=1.0/(1+exp(-probtemp))
            endif

            der=der+probtemp*(y(i,kv)-y(j,kv))*(x(i,lv)-x(j,lv))
            der2=der2+probtemp*(1-probtemp)*((y(i,kv)-y(j,kv))*(x(i,lv)-x(j,lv)))**2

          enddo
        enddo

        delta=(der+der2*theta(kv,lv))*2/(n-1)

        if(delta>=lam) then
          theta(kv,lv)=(delta-lam)/der2
          Sel(kv,lv)=1
        elseif(delta<=-lam) then
          theta(kv,lv)=(delta+lam)/der2
          Sel(kv,lv)=1
        else
          theta(kv,lv)=0
          Sel(kv,lv)=0
        endif

        if(Selold==Sel(kv,lv).and.abs(thetaold-theta(kv,lv))<eps) then
          !write(*,*)kv,lv,irep,abs(thetaold-theta(kv,lv))
          exit ITERATION
        endif

      enddo ITERATION

    enddo
  enddo

  !write(*,*)sum(Sel)

end subroutine penpwlkh
end module pwPENLKH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of pairwise likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
