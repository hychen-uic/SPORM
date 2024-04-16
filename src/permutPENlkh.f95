!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Network selection
!!       by penalized permutation likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------
!    Network selection
!--------------------------------------------------------------------------
! Input
!   dat(n,node)==> input data
!   ng ==> number of groups the nodes are divided.
!   group(ng) ==> number of nodes in each group
!   nlam ==> number of different penalty values
!   lambda(nlma) ==> different penalty values
! Output
!   Selnet==> selection outcomes Selnet(nlam,node,node)
!-------------------------------------------------------------------------
!######################################################################

    module permutPENLKH
    use, intrinsic :: iso_c_binding
    implicit none
    public :: NetworkSelectbyPM

    contains

!######################################################################

 subroutine NetworkSelectbyPM(dat,n,node,group,ng,lam,nlam,niter,eps,selnet,&
                              burnin,nsamp,sintv,maxcyc) bind(C, name = "networkselectbypm_")

   implicit none
   ! integer, parameter:: dp = selected_real_kind(15, 307)

   integer (C_INT) n,node,nlam,niter,ng !nlam is the number of lambda values to search
   integer (C_INT) burnin,nsamp,sintv
   real (C_DOUBLE) dat(n,node)
   integer (C_INT) group(ng),selnet(node,node,nlam)
   real (C_DOUBLE) lam(nlam)

   real (C_DOUBLE) eps,delta

   real (C_DOUBLE) theta(node,node),thetaold(node,node)
   !real(kind=dp),allocatable:: estv(:,:)

   integer (C_INT) maxcyc
   integer (C_INT) sample(nsamp,n,ng),tcyc(ng-1) ! sample of permutations
   real (C_DOUBLE) lowrate

   integer (C_INT) q,j,irep,i,k

   q=node*(node-1)/2 !Maximum number of parameters
   do j=1,ng
     q=q-group(j)*(group(j)-1)/2
   enddo

   do j=1, nlam
     selnet(:,:,j)=0
     theta=0
     tcyc=maxcyc
 !write(*,*)j,lam(j)

     IT:do irep=1,niter

 !write(*,*)irep,tcyc

       thetaold=theta
       ! 1.draw permutations for approximate permutation likelihood
       !call MSPwmgpen(dat,n,node,group,ng,selnet(:,:,j),theta,sample,nsamp,burnin,sintv,maxcyc)
       call MSPwmgpen(dat,n,node,group,ng,selnet(:,:,j),theta,sample,nsamp,burnin,sintv,maxcyc,tcyc,lowrate)
       maxcyc=maxval(tcyc) !adaptive steps sizes for sampling.
       ! 2. selection and estimation
       call SLAwmg(dat,n,node,group,ng,selnet(:,:,j),sample,nsamp,theta,q,lam(j),niter)
       ! 3. comparison for stop decision
       delta=sum(abs(theta-thetaold))
       if(delta<eps.and.irep>1) then
         exit IT
       endif

     enddo IT


     do k=1,node
       write(*,20)(selnet(i,k,j),i=1,node)
       20 Format(1x,50(I2,1x))
     enddo

   enddo

!   deallocate(estv)

 end subroutine NetworkSelectbyPM


!---------------------------------------------------
!  1.2. Laplace approximation for
!        a single group of outcomes
! Output a mean vector and a variance matrix
!---------------------------------------------------
subroutine SLAwmg(y,n,p,group,ng,network,sample,nsamp,theta,q,lambda,niter)
   !q=> total number of parameters to be estimated.
   !    q=p(p-1)/2-sumwk {group(k)(group(k)-1)/2}

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,p,nsamp,ng,q,niter
   integer group(ng),network(p,p),sample(nsamp,n,ng)
   real(kind=dp) y(n,p),lambda
   real(kind=dp) theta(p,p)!,estv(q,q)

   real(kind=dp) yy0(q),yy(nsamp,q),tmyy(q),myy(q),vyy(q,q)

   integer ind,vnet(q) ! vectorized network
   real(kind=dp) vtheta(q) ! vectorized block means

   integer i,j,k,ii,jj,kk,m,s,t

   m=p*(p-1)/2
   do j=1,ng
     m=m-group(j)*(group(j)-1)/2
   enddo
   if(m/=q) then
     !write(*,*) 'dimension does not match, stop excute LAwmg.'
     return
   endif

   ! vectorize theta to vtheta and compute statistics
    
   ind=0 !index non-zero positions
   kk=0
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
          ! vectorize theta
         vtheta(ind)=theta(s,t)!vectorize the non-zero theta
         vnet(ind)=network(s,t)
          ! compute statistic
         yy0(ind)=sum(y(:,s)*y(:,t))
         do ii=1,nsamp
           yy(ii,ind)=sum(y(sample(ii,:,k),s)*y(sample(ii,:,j),t))
         enddo
       enddo
       enddo
       jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

   do k=1,q ! calculate scores and second derivatives
     tmyy(k)=sum(yy(:,k))/nsamp
     do j=1,k
       vyy(k,j)=sum(yy(:,k)*yy(:,j))/nsamp-tmyy(k)*tmyy(j)
       vyy(j,k)=vyy(k,j)
     enddo
     myy(k)=yy0(k)-tmyy(k)
   enddo

   !call PDmatinv(vyy,q) !perform submatrix inversion

   !block-wise and column-wise conversion
   !call iUTvecF(network,p,group,ng,vnet,q) ! convert a block-wise matrix to vector
        ! network==>vnet
   Call QuadP(myy,-vyy,vtheta,vnet,q,niter,lambda) ! l_1 penalized selection

   !call iUTmatF(network,p,group,ng,vnet,q) ! convert a vector to bloack-wise matrix
        ! vnet==>network
   ind=0 !index non-zero positions
   kk=0
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
         network(s,t)=vnet(ind)
         network(t,s)=network(s,t)
       enddo
       enddo
       jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

end subroutine SLAwmg


!------------------------------------------------------------------------
!  1.3. Quadratic programming for estimating theta parameters
!  Input: mu (mean vector) and sig (the variance matrix),
!         lambda (penalty), and theta0.
!  Output: Parameter estimator for theta
!  The problem: min_{theta} mu*(theta-theta0)
!                 +0.5*(theta-theta0)*Sig*(theta-theta0)+lambda*|theta|
!  Given the current theta and current selected variables in SEL,
!       update both beta and SEL
!-------------------------------------------------------------------------
subroutine QuadP(mu,sig,gam,sel,nd,niter,lambda)
   implicit none
   integer nd,sel(nd),selold(nd),seld(nd),niter
   double precision mu(nd),sig(nd,nd),gam(nd),lambda

   integer k,irep,jv,nchange
   double precision delta,gam0(nd)

   gam0=gam ! record the initial value
   selold=sel  ! keep old selected set

   do irep=1,niter

     nchange=0  ! set change to 0.
     seld=sel  ! keep old selected set
     do jv=1,nd ! go through all the variables

       delta=mu(jv)-sig(jv,jv)*gam0(jv)
       do k=1,nd
         delta=delta+sig(jv,k)*(gam(k)*sel(k)-gam0(k)*selold(k))
       enddo
       delta=delta-sig(jv,jv)*(gam(jv)*sel(jv)-gam0(jv)*selold(jv)) ! remove the extra item should not be included

       if(delta>lambda) then
         gam(jv)=-(delta-lambda)/sig(jv,jv)
         Sel(jv)=1
       elseif(delta<-lambda) then
         gam(jv)=-(delta+lambda)/sig(jv,jv)
         Sel(jv)=1
       else
         gam(jv)=0
         Sel(jv)=0
       endif

     enddo

     nchange=sum(abs(sel-seld))
     if(nchange==0) then
       exit
     endif

   enddo

end subroutine QuadP

!-----------------------------------------------------------------------
 !  Metropolis sampling of permutations for multiple outcome groups
 !-----------------------------------------------------------------------
 subroutine MSPwmgpen(y,n,p,group,ng,network,theta,sample,nsamp,burnin,sintv,maxcyc,tcyc,lowrate)
    ! y==> y(n,p) outcomes. The last group of outcomes act as covariates
    ! ng==> the number of groups the p varaibles devided
    ! group==>group(np): the number of components in each group. Sum(group)=p
    ! theta==> theta[p,p] (symmetric) odds ratio parameters in the log-bilinear model
    ! n==> sample size of the observed data
    ! p==> dimension of covariates
    ! sample==> sample(nsamp,n,p): a sample of the permutations
    ! nsamp==> sample size of the generated permutation sample
    ! burnin==> length of burn in Markov chain
    ! sintv==> interval of sample selection after burn-in period

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,p,nsamp,burnin,sintv,maxcyc
   integer ng,group(ng),tcyc(ng-1),network(p,p)
   integer sample(nsamp,n,ng)! sampling permutations
   real(kind=dp) y(n,p)
   real(kind=dp) theta(p,p)

   integer isamp
   integer a(maxcyc),perm(n,ng)
   integer i,j,k,s,t,jj,kk,ll,itemp
   real(kind=dp) tot,runif,rate(ng-1),lowrate

!1. generate burnin sample and discard
   do i=1,n
     perm(i,1:ng)=i
   enddo ! initialize permutation

   do i=1,burnin !initial MCMC runs
     kk=0
     do k=1,ng-1 !Go through the first ng-1 groups of outcomes
                 !Permute within each group and decide if a permutation(sample)
                 !is drawn using the MCMC

       !Call SampleKfromN(a,maxcyc,n)
       Call SampleKfromN(a(1:tcyc(k)),tcyc(k),n)

       tot=0
       jj=0
       do j=1,ng
         if(j/=k) then
           do s=kk+1,kk+group(k)
           do t=jj+1,jj+group(j)
             if(network(s,t)/=0) then
               !do ll=1,maxcyc-1
               !  tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
               !          (y(perm(a(maxcyc),j),t)-y(perm(a(ll),j),t))*theta(s,t)
               !enddo
               do ll=1,tcyc(k)-1
                 tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                       (y(perm(a(tcyc(k)),j),t)-y(perm(a(ll),j),t))*theta(s,t)
               enddo
               !do ll=1,int(maxcyc/2) !for multiple cycles with length 2
               !  tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
               !          (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
               !enddo
             endif
           enddo
           enddo
         endif

         jj=jj+group(j)
       enddo
       kk=kk+group(k)


       Call random_number(runif)
       if(log(runif)<tot) then !Metropolis move
         itemp=perm(a(1),k) ! for one cycle with length greater than 1
         perm(a(1:(tcyc(k)-1)),k)=perm(a(2:tcyc(k)),k)
         perm(a(tcyc(k)),k)=itemp

         rate(k)=rate(k)+1
       endif

     enddo
   enddo

write(*,*)"rate=",rate*1.0/burnin
   do k=1,ng-1
     if((rate(k)*1.0/burnin)<lowrate.and.tcyc(k)>=2) then
       tcyc(k)=tcyc(k)-1
     endif
   enddo
write(*,*)"cycle=",tcyc

   isamp=0
   do while(isamp<nsamp) !MCMC runs for drawing sample
     do i=1,sintv ! sampling interval
       kk=0
       do k=1,ng-1 !permute the first ng-1 groups of outcomes

         Call SampleKfromN(a(1:tcyc(k)),tcyc(k),n)

         tot=0
         jj=0
         do j=1,ng
           if(j/=k) then
             do s=kk+1,kk+group(k)
             do t=jj+1,jj+group(j)
               if(network(s,t)/=0) then
                 do ll=1,tcyc(k)-1
                   tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                           (y(perm(a(tcyc(k)),j),t)-y(perm(a(ll),j),t))*theta(s,t)
                 enddo
                 !do ll=1,int(maxcyc/2) !for multiple cycles with length 2
                 !  tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
                 !          (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
                 !enddo
               endif
             enddo
             enddo
           endif

           jj=jj+group(j)
         enddo
         kk=kk+group(k)

         call random_number(runif)
         if(log(runif)<tot) then !Metropolis move
           itemp=perm(a(1),k)
           perm(a(1:(tcyc(k)-1)),k)=perm(a(2:tcyc(k)),k)
           perm(a(tcyc(k)),k)=itemp

            ! do ll=1,int(maxcyc/2)! for multiple cycles with length 2
            !   itemp=perm(a(2*ll-1),k)
            !   perm(a(2*ll-1),k)=perm(a(2*ll),k)
            !   perm(a(2*ll),k)=itemp
            ! enddo
         endif

       enddo
     enddo

     isamp=isamp+1
     sample(isamp,:,:)=perm
   enddo

   !write(*,*)"Out of ",sintv*nsamp, " attempts ", nmove, " moved."

 end subroutine MSPwmgpen


 !-----------------------------------------------------------------------
 !  Metropolis sampling of permutations for multiple outcome groups
 !-----------------------------------------------------------------------
 subroutine MSPwmgpen0(y,n,p,group,ng,network,theta,sample,nsamp,burnin,sintv,maxcyc)
    ! y==> y(n,p) outcomes. The last group of outcomes act as covariates
    ! ng==> the number of groups the p varaibles devided
    ! group==>group(np): the number of components in each group. Sum(group)=p
    ! theta==> theta[p,p] (symmetric) odds ratio parameters in the log-bilinear model
    ! n==> sample size of the observed data
    ! p==> dimension of covariates
    ! sample==> sample(nsamp,n,p): a sample of the permutations
    ! nsamp==> sample size of the generated permutation sample
    ! burnin==> length of burn in Markov chain
    ! sintv==> interval of sample selection after burn-in period

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,p,nsamp,burnin,sintv,maxcyc
   integer ng,group(ng),network(p,p)
   integer sample(nsamp,n,ng)! sampling permutations
   real(kind=dp) y(n,p)
   real(kind=dp) theta(p,p)

   integer isamp
   integer a(maxcyc),perm(n,ng)
   integer i,j,k,s,t,jj,kk,ll,itemp
   real(kind=dp) tot,runif

!1. generate burnin sample and discard
   do i=1,n
     perm(i,1:ng)=i
   enddo ! initialize permutation

   do i=1,burnin !initial MCMC runs
     kk=0
     do k=1,ng-1 !Go through the first ng-1 groups of outcomes
                 !Permute within each group and decide if a permutation(sample)
                 !is drawn using the MCMC

       Call SampleKfromN(a,maxcyc,n)
       tot=0
       jj=0
       do j=1,ng
         if(j/=k) then
           do s=kk+1,kk+group(k)
           do t=jj+1,jj+group(j)
             if(network(s,t)/=0) then
               !do ll=1,maxcyc-1
               !  tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
               !          (y(perm(a(maxcyc),j),t)-y(perm(a(ll),j),t))*theta(s,t)
               !enddo
               do ll=1,int(maxcyc/2) !for multiple cycles with length 2
                 tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
                         (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
               enddo
             endif
           enddo
           enddo
         endif

         jj=jj+group(j)
       enddo
       kk=kk+group(k)

       Call random_number(runif)
       if(log(runif)<tot) then !Metropolis move
         !itemp=perm(a(1),k) ! for one cycle with length greater than 1
         !perm(a(1:(maxcyc-1)),k)=perm(a(2:maxcyc),k)
         !perm(a(maxcyc),k)=itemp
         do ll=1,int(maxcyc/2)! for multiple cycles with length 2
           itemp=perm(a(2*ll-1),k)
           perm(a(2*ll-1),k)=perm(a(2*ll),k)
           perm(a(2*ll),k)=itemp
         enddo
       endif

     enddo
   enddo

   !write(*,*)"Out of ",burnin, " attempts ", nmove, " moved."


   isamp=0
   do while(isamp<nsamp) !MCMC runs for drawing sample
     do i=1,sintv ! sampling interval
       kk=0
       do k=1,ng-1 !permute the first ng-1 groups of outcomes

         Call SampleKfromN(a,maxcyc,n)
         tot=0
         jj=0
         do j=1,ng
           if(j/=k) then
             do s=kk+1,kk+group(k)
             do t=jj+1,jj+group(j)
               if(network(s,t)/=0) then
                 !do ll=1,maxcyc-1
                 !  tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                 !          (y(perm(a(maxcyc),j),t)-y(perm(a(ll),j),t))*theta(s,t)
                 !enddo
                 do ll=1,int(maxcyc/2) !for multiple cycles with length 2
                   tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
                           (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
                 enddo
               endif
             enddo
             enddo
           endif

           jj=jj+group(j)
         enddo
         kk=kk+group(k)

         call random_number(runif)
         if(log(runif)<tot) then !Metropolis move
           !itemp=perm(a(1),k)
           !perm(a(1:(maxcyc-1)),k)=perm(a(2:maxcyc),k)
           !perm(a(maxcyc),k)=itemp
           do ll=1,int(maxcyc/2)! for multiple cycles with length 2
             itemp=perm(a(2*ll-1),k)
             perm(a(2*ll-1),k)=perm(a(2*ll),k)
             perm(a(2*ll),k)=itemp
           enddo
         endif

       enddo
     enddo

     isamp=isamp+1
     sample(isamp,:,:)=perm
   enddo

   !write(*,*)"Out of ",sintv*nsamp, " attempts ", nmove, " moved."

 end subroutine MSPwmgpen0

!-----------------------------------------------------------
! 2.2. Sampling K indices from {1,...,n} without replacement
!----------------------------------------------------------
 subroutine SampleKfromN(sample,K,N)

   implicit none

   integer K,N
   integer Sample(K) ! the sampled indices

   integer i,j,index
   integer CSO(K) ! Current ordered sample: ordered samples that have already generated
   double precision runif(K)

   Call random_number(runif)
   sample(1)=int(n*runif(1))+1
   CSO(1)=sample(1)
   do i=2,K
     sample(i)=int((n-i+1)*runif(i))+1
     index=1
     do j=1,i-1 ! find the right place for the current sample
       if(sample(i)>=CSO(j)) then
         sample(i)=sample(i)+1
         index=index+1
       endif
     enddo
     if(index==i) then !adjust CSO
       CSO(i)=sample(i)
     else
       do j=i,index+1,-1
         CSO(j)=CSO(j-1)
       enddo
       CSO(index)=sample(i)
     endif
   enddo
   !sample=CSO ! this is needed if sample is ordered.

 end subroutine SampleKfromN

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
! Simple inverse by Cholestry decomposition
!---------------------------------------------
subroutine PDmatinv(a,n)
      integer n
      DOUBLE PRECISION a(n,n)

      integer i,j,k
      double precision tot,b(n)

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
                 ! if(tot<=0) then
                 !    write(*,*)j,j,'  choldc failed'
                 ! endif
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

 !--------------------------------------------------
!  Convert network structure (block-wise) into
!     a vector to output
!--------------------------------------------------
 subroutine iUTvecF(network,p,group,ng,vnet,q)
    implicit none

    integer p,ng,q
    integer network(p,p),group(ng),vnet(q)

    integer j,k,jj,kk,count,s,m

    m=p*(p-1)/2 !Maximum number of parameters
    do j=1,ng
      m=m-group(j)*(group(j)-1)/2
    enddo
    if(m/=q) then
      write(*,*)"Dimension q does not match with the dimensions of the block matrices"
      return
    endif

    count=0
    kk=0
    do k=1,ng-1
      jj=kk+group(k)
      do j=(k+1),ng
        do s=jj+1,jj+group(j)
          vnet(count+1:count+group(k))=network(kk+1:kk+group(k),s)
          count=count+group(k)
        enddo
        jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

end subroutine iUTvecF

!--------------------------------------------------
!  Convert the vector structure into
!     a block-wise matrix (inverse of iUTvecF)
!--------------------------------------------------
 subroutine iUTmatF(network,p,group,ng,vnet,q)
    implicit none

    integer p,ng,q
    integer network(p,p),group(ng),vnet(q)

    integer j,k,jj,kk,count,s,m!,t

    m=p*(p-1)/2 !Maximum number of parameters
    do j=1,ng
      m=m-group(j)*(group(j)-1)/2
    enddo
    if(m/=q) then
      write(*,*)"Dimension q does not match with the dimensions of the block matrices"
      return
    endif

    count=0
    kk=0
    do k=1,ng-1
      jj=kk+group(k)
      do j=(k+1),ng
        do s=jj+1,jj+group(j)
          network(kk+1:kk+group(k),s)=vnet(count+1:count+group(k))
          network(s,kk+1:kk+group(k))=network(kk+1:kk+group(k),s)
          count=count+group(k)
        enddo
        jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

end subroutine iUTmatF
end module permutPENLKH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    End of network selection by
!!        penalized permutation likelihood approach
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
