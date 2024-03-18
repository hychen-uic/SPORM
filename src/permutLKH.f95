!-----------------------------------------------------
! Analysis by SORM with grouped multiple outcomes
!-----------------------------------------------------
!######################################################################

    module permutLKH
    use, intrinsic :: iso_c_binding
    implicit none
    public :: Analysiswmg

    contains

!######################################################################
subroutine Analysiswmg(dat,n,p,group,ng,vtheta,estv,q,eps,converge,niter,nlag,burnin,nsamp, &
                       sintv,maxcyc,theta2,steplen,nstep) bind(C, name = "analysiswmg_")
 ! dat(n,p)=> outcomes and covariates
 ! ng==> the number of groups the p varaibles devided
 ! group==>group(np): the number of components in each group. Sum(group)=p
 ! niter=>Total number of iterations
 ! nlag=> number of iteration to be averaged for convergence assessment
 ! burnin=> Gibbs sampler burnin runs
 ! nsamp=> sample size of the draws
 ! sintv=> interval between two sample draws
 ! steplen=> control the step length in the parameter update
 !             steplen=1<=> no control; steplen<1 <=> shorten step size
 ! nstep=> control the number of steps before a change of the step size occurs 

!   parameter(n=200,p=5,nstart=31,niter=50,allrep=100)
!   parameter(burnin=500000,Nsamp=2000,Sintv=10000)

   implicit none
   !integer, parameter:: dp = selected_real_kind(15, 307)

   integer (C_INT) n,p,ng,q,converge,niter,nlag,burnin,nsamp,sintv,maxcyc,nstep
   integer (C_INT) group(ng),tcyc(ng-1)

   real(C_DOUBLE) dat(n,p)
   real(C_DOUBLE) theta(p,p),theta2(niter,p,p)
   real(C_DOUBLE) vtheta(q),estv(q,q),estv2(niter,q,q)  !q=p*(p-1)-sum_k group(k)(group(k)-1)/2
   real(C_DOUBLE) mtheta,vartheta,vmax,eps,lowrate,steplen(2),stepsize

   integer (C_INT) i,j,k,irep,kk,jj,s,t,last,ind,ESS
   integer (C_INT) sample(nsamp,n,ng)

   do i=1,n  !initialize
      sample(1:nsamp,i,1:ng)=i
   enddo

   ind=0 !matrize vtheta to matrix theta
   kk=0
   theta=0 !initialize theta, not all components of theta will be updated
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
         theta(s,t)=vtheta(ind) !matrize vtheta to theta
       enddo
       enddo
       jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

   !open(unit=5,file="result5",status="unknown")

   estv=0
   tcyc=maxcyc
   IT:do irep=1,niter

     !write(*,*)irep,niter
     !lowrate=ESS*1.0/(sintv*nsamp) !at least to have effective sample size 2000.
     if(irep<100) then
       lowrate=1e-5
     elseif(irep<200) then
       lowrate=1e-4
     else
       lowrate=1e-3
     endif

     Call MSPwmg(dat,n,p,group,ng,theta,sample,nsamp,burnin,sintv,maxcyc,tcyc,lowrate)
     maxcyc=maxval(tcyc) !adaptive steps sizes for sampling.

     if(irep<nstep) then
       stepsize=steplen(1)
     else
       stepsize=steplen(2)
     endif
     Call LAwmg(dat,n,p,group,ng,sample,nsamp,theta,estv,q,stepsize)

    ! write(5,10)irep,((theta(k,j),k=1,j-1),j=2,p)
    ! 10 format(1x,I5, 2x, 100(f12.4,2x))

     theta2(irep,:,:)=theta
     estv2(irep,:,:)=estv

     if(irep>=(nlag+20)) then ! check convergence
       vmax=0
       kk=0
       do k=1,ng-1
         jj=kk+group(k)
         do j=(k+1),ng
           do s=kk+1,kk+group(k)
             do t=jj+1,jj+group(j)
               mtheta=sum(theta2((irep-nlag+1):irep,s,t))/nlag
               vartheta=sum(theta2((irep-nlag+1):irep,s,t)**2)/nlag
               vartheta=vartheta-mtheta**2
               if(vartheta/(1+mtheta*mtheta)>vmax) then
                 vmax=vartheta/(1+mtheta*mtheta)
               endif
             enddo
           enddo
           jj=jj+group(j)
         enddo
         kk=kk+group(k)
       enddo

       if(sqrt(vmax)<eps) then
         converge=1
         last=irep
         exit IT
       elseif(irep==niter)then
         converge=0
         last=niter
       endif
     endif

   enddo IT

   !note irep is the last iteration index which should be used when loop is interrupted
   theta=0
   do j=1,p-1
     do k=j+1,p
       theta(j,k)=sum(theta2((last-nlag+1):last,j,k))/nlag
       theta(k,j)=theta(j,k)
     enddo
   enddo

   estv=0
   do j=1,q
     do k=1,j
       estv(j,k)=sum(estv2((last-nlag+1):last,j,k))/nlag
       estv(k,j)=estv(j,k)
     enddo
   enddo

   ind=0 !convert matrix theta back into a vector
   kk=0
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
         vtheta(ind)=theta(s,t) !vectorize theta
       enddo
       enddo
       jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

end subroutine Analysiswmg

!-------------------------------------------------------------
! 1.4. Laplace approximation for multiple groups of outcomes
!-------------------------------------------------------------
subroutine LAwmg(y,n,p,group,ng,sample,nsamp,theta,estv,q,stepsize)
   !q=> total number of parameters to be estimated.
   !    q=p(p-1)/2-sumwk {group(k)(group(k)-1)/2}

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer n,p,nsamp,ng,group(ng),q
   real(kind=dp) y(n,p)
   real(kind=dp) theta(p,p),estv(q,q)
   integer sample(nsamp,n,ng) ! sampled permutations of y groups

   real(kind=dp) yy0(q),yy(nsamp,q),tmyy(q),myy(q),vyy(q,q)

   integer ind ! vectorized network
   real(kind=dp) vtheta(q) ! vectorized block means

   integer i,j,k,ii,jj,kk,m,s,t
   real(kind=dp) stepsize

   m=p*(p-1)/2
   do j=1,ng
     m=m-group(j)*(group(j)-1)/2
   enddo
 !  if(m/=q) then
 !    !write(*,*) 'dimension does not match, stop excute LAwmg.'
 !    return
 !  endif

   ! vectorize theta to vtheta and compute statistics

   ind=0 !index non-zero positions
   kk=0
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
           !vectorize theta
         vtheta(ind)=theta(s,t)
           !compute statiatic
         yy0(ind)=sum(y(:,s)*y(:,t)) ! for a given parameter, sum over all subjects
         do ii=1,nsamp
           yy(ii,ind)=sum(y(sample(ii,:,k),s)*y(sample(ii,:,j),t)) !sum over all subjects
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

   call PDmatinv(vyy,q) !perform submatrix inversion

   vtheta=vtheta+stepsize*matmul(vyy,myy) !one-step update

 ! Laplace Apprx estimated variance is Sigma^(-1)exp(-theta*yxmean/2)

   estv=vyy

   ind=0 !convert vtheta back to matrix position
   kk=0
   theta=0 !initialize theta, not all components of theta will be updated
   do k=1,ng-1
     jj=kk+group(k)
     do j=(k+1),ng
       do t=jj+1,jj+group(j)
       do s=kk+1,kk+group(k) !t,s order determined the vectorization
         ind=ind+1
         theta(s,t)=vtheta(ind) !vectorize the non-zero theta
         theta(t,s)=theta(s,t)
       enddo
       enddo
       jj=jj+group(j)
     enddo
     kk=kk+group(k)
   enddo

end subroutine LAwmg

 !-----------------------------------------------------------------------
 !  Metropolis sampling of permutations for multiple outcome groups
 !-----------------------------------------------------------------------
 subroutine MSPwmg(y,n,p,group,ng,theta,sample,nsamp,burnin,sintv,maxcyc,tcyc,lowrate)
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
   integer ng,group(ng)
   integer sample(nsamp,n,ng)! sampling permutations
   real(kind=dp) y(n,p)
   real(kind=dp) theta(p,p)

   integer isamp
   integer a(maxcyc),perm(n,ng)
   integer i,j,k,s,t,jj,kk,ll,itemp
   real(kind=dp) tot,runif,lowrate

   integer rate(ng-1),tcyc(ng-1)

!1. generate burnin sample and discard
   do k=1,ng
     do i=1,n
       perm(i,k)=i
     enddo
   enddo ! initialize permutation

   rate=0
   do i=1,burnin
     kk=0
     do k=1,ng-1 !permute the first ng-1 groups of outcomes

       !Call SampleKfromN(a,maxcyc,n)
       Call SampleKfromN(a(1:tcyc(k)),tcyc(k),n)

       tot=0
       jj=0
       do j=1,ng
         if(j/=k) then
           do s=kk+1,kk+group(k)
           do t=jj+1,jj+group(j)
             do ll=1,tcyc(k)-1
               !tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
               !        (y(perm(a(maxcyc),j),t)-y(perm(a(ll),j),t))*theta(s,t)
               tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                       (y(perm(a(tcyc(k)),j),t)-y(perm(a(ll),j),t))*theta(s,t)
             enddo
             !do ll=1,int(maxcyc/2) !for multiple cycles with length 2
             !  tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
             !          (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
             !enddo

           enddo
           enddo
         endif!

         jj=jj+group(j)
       enddo
       kk=kk+group(k)

       Call random_number(runif)
       if(log(runif)<tot) then !Metropolis move
         itemp=perm(a(1),k) ! for one cycle with length greater than 1
         !perm(a(1:(maxcyc-1)),k)=perm(a(2:maxcyc),k)
         !perm(a(maxcyc),k)=itemp
         perm(a(1:(tcyc(k)-1)),k)=perm(a(2:tcyc(k)),k)
         perm(a(tcyc(k)),k)=itemp

         rate(k)=rate(k)+1

         !do ll=1,int(maxcyc/2)! for multiple cycles with length 2
         !  itemp=perm(a(2*ll-1),k)
         !  perm(a(2*ll-1),k)=perm(a(2*ll),k)
         !  perm(a(2*ll),k)=itemp
         !enddo
       endif

     enddo
   enddo

   !write(*,*)"rate=",rate*1.0/burnin
   do k=1,ng-1
     if((rate(k)*1.0/burnin)<lowrate.and.tcyc(k)>=2) then
       tcyc(k)=tcyc(k)-1
     endif
   enddo
   !write(*,*)"cycle=",tcyc

   isamp=0
   do while(isamp<nsamp)
     do i=1,sintv ! sampling interval
       kk=0
       do k=1,ng-1 !permute the first ng-1 groups of outcomes

         !Call SampleKfromN(a,maxcyc,n)
         Call SampleKfromN(a(1:tcyc(k)),tcyc(k),n)

         tot=0 ! compute move probability ratio
         jj=0
         do j=1,ng
           if(j/=k) then
             do s=kk+1,kk+group(k)
             do t=jj+1,jj+group(j)
               do ll=1,tcyc(k)-1
                 !tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                 !        (y(perm(a(maxcyc),j),t)-y(perm(a(ll),j),t))*theta(s,t)
                 tot=tot+(y(perm(a(ll),k),s)-y(perm(a(ll+1),k),s))* &
                         (y(perm(a(tcyc(k)),j),t)-y(perm(a(ll),j),t))*theta(s,t)
               enddo
               !do ll=1,int(maxcyc/2) !for multiple cycles with length 2
               !  tot=tot+(y(perm(a(2*ll-1),k),s)-y(perm(a(2*ll),k),s))* &
               !          (y(perm(a(2*ll),j),t)-y(perm(a(2*ll-1),j),t))*theta(s,t)
               !enddo
             enddo
             enddo
           endif!

           jj=jj+group(j)
         enddo
         kk=kk+group(k)

         call random_number(runif)
         if(log(runif)<tot) then !Metropolis move
           itemp=perm(a(1),k)
           !perm(a(1:(maxcyc-1)),k)=perm(a(2:maxcyc),k)
           !perm(a(maxcyc),k)=itemp
           perm(a(1:(tcyc(k)-1)),k)=perm(a(2:tcyc(k)),k)
           perm(a(tcyc(k)),k)=itemp
           !do ll=1,int(maxcyc/2)! for multiple cycles with length 2
           !  itemp=perm(a(2*ll-1),k)
           !  perm(a(2*ll-1),k)=perm(a(2*ll),k)
           !  perm(a(2*ll),k)=itemp
           !enddo
         endif

       enddo
     enddo

     isamp=isamp+1
     sample(isamp,:,:)=perm
   enddo

   !write(*,*)"Out of ",sintv*nsamp, " attempts ", nmove, " moved."

 end subroutine MSPwmg

!-----------------------------------------------------------
! sampling K indices from {1,...,n} without replacement
!----------------------------------------------------------
 subroutine SampleKfromN(sample,K,N)

   implicit none
   integer, parameter:: dp = selected_real_kind(15, 307)

   integer K,N
   integer Sample(K) ! the sampled indices

   integer i,j,index
   integer CSO(K) ! Current ordered sample: ordered samples that have already generated
   real(kind=dp) runif(K)

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

!--------------------------------------------
! Simple inverse by Cholestry decomposition
!---------------------------------------------
subroutine PDmatinv(a,n)

      implicit none
      integer, parameter:: dp = selected_real_kind(15, 307)

      integer n
      real(kind=dp) a(n,n)

      integer i,j,k
      real(kind=dp) tot,b(n)

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
end module permutLKH
