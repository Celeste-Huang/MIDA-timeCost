!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

module MCMC
    use DALEC_deciduous
    implicit none
    
    character(len=150)::obsFileName='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\obsNEE.txt'
    character(len=150)::paramFileName='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\param.txt'
    integer:: paramNum=21
    integer:: time=5844
    integer:: D=50
    real(kind=8),dimension(21)::cBest,c,cmin,cmax,cNew
    integer:: nsimu=20000
    real(kind=8),dimension(21,20001):: c_record
    real(kind=8),dimension(20001):: J_record
    integer:: record
    real(kind=8),dimension(5844):: simuBest,obs,simu
    real(kind=8):: var,J_new,J_default=20000000
    character(len=150)::outJ='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\output\mismatch_accepted.txt'
    character(len=150)::outC='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\output\param_accepted.txt'
    character(len=150)::outCBest='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\output\param_best.txt'
    character(len=150)::outSimuBest='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\output\simu_best.txt'
    character(len=150)::outRecord='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\output\record.txt'
    contains
    
  subroutine readObsData()
      real(kind=8) :: mean,vartemp
      integer :: i
      open(100,file=obsFileName)
      read(100,*) (obs(i), i=1,time)
      close(100)
      mean=SUM(obs)/time
      vartemp=(obs(1)-mean)**2
      do i=2,time
          vartemp=vartemp+(obs(i)-mean)**2
      end do
      var=vartemp/(time-1)
  end subroutine readObsData

  subroutine readParam()
      integer :: i
      open(100,file=paramFileName)
      read(100,*) (cmin(i), i=1,paramNum)
      read(100,*) (cmax(i), i=1,paramNum)
      read(100,*) (c(i), i=1,paramNum)
      close(100)
  end subroutine readParam
  
  logical function is_qualified(cNew)
      real(kind=8),dimension(:),intent(in) :: cNew
      integer :: i
      is_qualified= .True.
      do i=1,paramNum
          if(.not. ((cmin(i)<cNew(i)) .and. (cmax(i)>cNew(i))))then
              print*, i,'th parameter fails to get a reasonable value'
          end if
          is_qualified = is_qualified .and. (cmin(i)<cNew(i)) .and. (cmax(i)>cNew(i))
      end do
  end function is_qualified

      subroutine genParam(c_opt,c_new)
          real(kind=8),dimension(21)::  c_opt
          real(kind=8),dimension(21):: c_new
          real(kind=8),dimension(21):: randn,randn1
          integer :: i
          do while(.true.)
             call random_number(randn1)
             randn1=randn1-0.5
             do i=1, paramNum
                c_new(i)=c_opt(i)+randn1(i)*(cmax(i)-cmin(i))/D
             enddo
             if (is_qualified(c_new)) exit 
          enddo
          !open(100,file='D:\Dropbox\00Work\MIDA-submit-GMD\review\Trade-offs\TradionalDA\paramValue1.txt')
          !write(100,'(F30.15)') (c_new(i), i=1,21)
          !close(100)
      end subroutine genParam  
      
  subroutine get_error(err)
      real(kind=8) :: err
      integer :: i
      err=0
      do i=1,time
          err=err+(obs(i)-simu(i))**2
      end do
      err=err/(2*var)
  end subroutine get_error
  
  subroutine mcmc_alg()
          real(kind=8),dimension(paramNum) :: randn
          real(kind=8) :: scaleVar, randn1,randn2,J_new,delta_J
          integer :: isimu,i,j
          record=1
          call random_number(randn1)

          c_record(:,1)=c
          J_record(1)=J_default
          do isimu=1,nsimu
              !! proposing step
              call genParam(c_record(:,record),cNew)
              call runModel(cNew,simu)
              call get_error(J_new)
              !print *,J_new
              !! moving step
              delta_J=J_record(record)-J_new
              call random_number(randn2)
              if(min(1.0,exp(delta_J))>randn2) then
                  record=record+1
                  c_record(:,record)=cNew
                  J_record(record)=J_new
              end if
              write(*, '("isimu = ", I5,"   upgraded = " , I5    "   mismatch = " , F30.15    "   acceptRatio = " , F30.15)') &
              &isimu, (record-1), J_new,(record-1.0)/isimu
          end do
      end subroutine mcmc_alg
      
      subroutine getBest()
          integer :: index(1)
          index=minloc(J_record(1:record))
          cBest=c_record(:,index(1))
          call runModel(cBest,simuBest)
      end subroutine getBest
     
      subroutine write_io_file()
        integer :: i,j,k
        character(len=1) :: M
        write (M,'(I1)') paramNum
        open(unit=100    ,file=trim(adjustl(outJ)))
        open(unit=101    ,file=trim(adjustl(outC)))
        open(unit=102    ,file=trim(adjustl(outRecord)))
        open(unit=103    ,file=trim(adjustl(outSimuBest)))
        open(unit=104    ,file=trim(adjustl(outCBest)))
        do i=101,record
            write(100,'(F30.15)') J_record(i)            
            write(101,'('//M//'(F30.15,'//'" "))') (c_record(j,i), j=1,paramNum)
        end do
        write(102,*) record-1
        write(103,*) (simuBest(j),j=1,time)
        write(104,*) cBest
        close(100)
        close(101)
        close(102)
        close(103)
        close(104)
    end subroutine write_io_file
     end module MCMC