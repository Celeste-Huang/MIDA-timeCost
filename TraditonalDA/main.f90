!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!
program main
    !use DALEC_deciduous
    use MCMC    
    implicit none 
    
    character(len=10)   :: time1,time2
    character(len=8)    :: thedate
    real(kind=4)    :: t1,t2
    call date_and_time(thedate,time1)
    write(*,*) 'date: ',thedate
    write(*,*) 'start time:',time1(1:2) // ':' // time1(3:4) // ':' // time1(5:10)
    call readObsData()
    call readParam()
    call mcmc_alg()   
    call getBest()
    call write_io_file()
    
    call date_and_time(thedate,time2)
    read(time1,*) t1
    read(time2,*) t2
    write(*,*) 'end time:',time2(1:2) // ':' // time2(3:4) // ':' // time2(5:10)
    write(*,*) 'time cost in Traditional DA:', t2-t1
    !call runDALEC()
    !call runModel()
end program main
    
