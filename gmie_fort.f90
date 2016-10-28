module gMIE_module
  implicit none
contains

  subroutine gMIE_core(wl_0,r_p,eper_m,mper_m,eper_p,mper_p,nang,Qsca,Qabs,Qext,Qext_opt,S1,S2)
  !=================================== subroutine gMIE_core =======================================
  ! Computes Mie coefficients a,b,c,d and far-field quantities Qsca,Qabs,Qext,Qext_opt,S1,S2
  ! for of an isotropic homogeneous sphere based on the MIE theory generalized for magnetic particle and abosorbing media.
  ! permittivity and permiability of particle (ε1  and μ1) and those of medium (ε and μ) are all complex scalar.

  !-----References-----
  !BH83: Bohren and Huffman 1983, Absorption and Scatteing of Light by Small Particles, John Wiley & Sons. Inc.
  !KER83: Kerker et al. 1983, Electromagnetic scattering by magnetic spheres, J.Opt.Soc.Am., vol. 73, 765-767.
  !MISH07: Mischchenko 2007, Electromagnetic scattering by a fixed finite object embedded in an absorbing medium, Optics Express,vol 15, 13188-13202.
  !--------------------

  !---- History ----
  !version 1.2: Nov. 23, 2015, by Nobuhiro Moteki

  !--- Notes ---
  ! refractive index of particle m_p is defined as sqrt(eper_p*mper_p)
  ! refractive index of medium m_m is defined as sqrt(eper_m*mper_m)
  ! relative refractive index of particle defined as m_p/m_m
  ! incident field is plane wave propagating into +z direction
  ! In absorbing media, far-field quantities is only Qext_opt. In nonabsorbing media, Qext = Qext_opt. See MISH07 for detail.
  ! Near-field quantities can be computed using the Mie coefficients a,b,c,d.


    !---INPUT ARGUMENTS---
    double precision wl_0 ! wavelength in vacuum (=c/w)
    double precision  r_p ! particle radius
    complex(kind(0d0)) eper_m ! complex permittivity of medium ε
    complex(kind(0d0)) mper_m ! complex permiability of medium μ
    complex(kind(0d0)) eper_p ! complex permittivity of particle ε1
    complex(kind(0d0)) mper_p ! complex permiability of particle μ1
    integer nang ! number of scattering angle between 0-180 deg for scattering matrix element S1 S2
    !---------------------

    !---OUTPUT ARGUMENTS(far-field quantities)---
    double precision Qsca !  scattering efficiency := Csca/(pi*r_p**2)
    double precision Qext ! extinction efficiency := Cext/(pi*r_p**2)
    double precision Qabs ! absorption efficiency := Cabs/(pi*r_p**2)
    double precision Qext_opt ! extinction efficiency based on optical theorem := Cext/(pi*r_p**2) defined as MISH07
    complex(kind(0d0)) S1(nang) ! (2,2) element of the 2*2 amplitude scattering matrix defined as BH83, Eq.3.12
    complex(kind(0d0)) S2(nang) ! (1,1) element of the 2*2 amplitude scattering matrix defined as BH83, Eq.3.12
    !---------------------

    !--- internal variables ---
    complex(kind(0d0)) m_m ! refractive index of medium
    complex(kind(0d0)) m_p ! complex refractive index of particle
    complex(kind(0d0)) m_r ! relative complex refractive index of particle (=m_p/m_m)
    double precision k0 ! free space wavenumber
    complex(kind(0d0)) k ! wavenumber in medium
    complex(kind(0d0)) x ! size parameter of particle with respect to the medium (=2*pi*r_p*m_m/wl_0)
    complex(kind(0d0)) y ! auxiliary parameter for numerical computation
    integer nstop ! number of Mie terms
    integer nmx ! number of terms for downward recurrence
    complex(kind(0d0)),allocatable :: DD(:) ! logarithmic derivative  (BH83 Eq.4.89)
    complex(kind(0d0)),allocatable :: psi(:),chi(:),qsi(:),psim(:),RR(:) !Reccati-Bessel functions (0:nstop)
    complex(kind(0d0)),allocatable :: a(:),b(:),c(:),d(:) !Mie coefficients (1:nstop)
    double precision ,allocatable :: fn1(:),fn2(:)
    double precision theta,dtheta,mu
    double precision, allocatable :: pie(:) ! angular function pie (BH83)
    double precision, allocatable :: tau(:) ! angular function tau (BH83)
    integer n,j ! loop index
    !--------------------------

    complex(kind(0d0)), parameter :: i=(0.d0,1.d0) ! imaginary number
    double precision, parameter :: pi=acos(-1.0d0)
    !============================================================================================

    k0=2*pi/wl_0 ! free space wavenumber
    m_m=sqrt(eper_m*mper_m) ! refractive index of medium
    m_p=sqrt(eper_p*mper_p) ! refractive index of particle
    k=m_m*k0 ! wavenumber in medium
    x=k*r_p ! size parameter of particle in medium
    m_r=m_p/m_m ! relative refractive index of particle
    y=x*m_r ! auxiliary parameter for numerical computation

    nstop=floor(abs(x)+4*abs(x)**0.3333+2) ! number of expansion terms for partial wave coefficients (BH83)
    allocate(chi(0:nstop),qsi(0:nstop))
    allocate(a(1:nstop),b(1:nstop),c(1:nstop),d(1:nstop))
    allocate(fn1(1:nstop),fn2(1:nstop))
    allocate(pie(1:nstop),tau(1:nstop))

    nmx=max(nstop,int(abs(y)))+15 ! number of terms for downward recurence calculations
    allocate(DD(1:nmx),psi(0:nmx),psim(0:nmx),RR(0:nmx))

    !---------------Logarithmic derivative DD calculated by downward recurrence-------------------
    ! beginning with initial value (0.,0.) at nmx
    ! y:=x*m_r argument of DD
    DD=0.d0
    do n=nmx,2,-1
        DD(n-1)=dble(n)/y-1.d0/(DD(n)+dble(n)/y)
    enddo
    DD(1:nstop)=reshape(DD,(/nstop/)) ! trim redundant elements
    !---------------------------------------------------------------------------------------------

    !-----------Reccati-Bessel function PSI(0:nstop) calculated by downward recurrence----------
    ! Reference: Mischenko et al. 2002, Scattering, Absorption and Emission of Light by Small Particles 3rd, pp.167-169
    ! x:=k*r_p argument of PSI
    RR=0.0d0; ! R(n):=PSI(n)/PSI(n-1)
    RR(nmx)=x/(2.d0*nmx+1.d0) ! starting value of downward recurrence
    do n=nmx-1,0,-1
        RR(n)=1.d0/((2.d0*n+1.d0)/x-RR(n+1)) ! R(n) := Rn
    enddo

    psi(0)=RR(0)*cos(x)
    do n=1,nstop
        psi(n)=RR(n)*psi(n-1) ! PSI(n) := PSIn
    enddo
    psi(0:nstop)=reshape(psi,(/nstop+1/)) !trim redundant elements
    !---------------------------------------------------------------------------------------------

    ! -------Reccati-Bessel function chi(0:nstop) calculated by upward recurrence------
    ! beginning with initial values chi(-1)=sin(x), chi(0)=-cos(x)
    ! Reference: Bphren and Huffman 1983 (Appendix A, p.478)
    ! CHI(x):= x*y(x) where y(x) is spherical bessel function of second kind
    ! This is contrast to the definition of BH83: CHI(x):= -x*y(x)
    ! x:=k*r_p argument of CHI
    chi(0)=-cos(x)
    chi(1)=(1.d0/x)*chi(0)-sin(x)
    do n=2,nstop
        chi(n)=((2.0d0*n-1.0d0)/x)*chi(n-1)-chi(n-2)
    enddo


    !---------------------------------------------------------------------------------------------

    do n=0,nstop
    qsi(n)=psi(n)+i*chi(n) ! Reccati-Bessel function of third kind := x*(j(x)+iy(x))
    enddo



    !-----------Reccati-Bessel function PSIM(0:nstop) calculated by downward recurrence----------
    ! Reference: Mischenko et al. 2002, Scattering, Absorption and Emission of Light by Small Particles 3rd, pp.167-169
    ! y:=x*m_r argument of PSI
    RR=0.0d0; ! R(n):=PSI(n)/PSI(n-1)
    RR(nmx)=y/(2.d0*nmx+1.d0) ! starting value of downward recurrence
    do n=nmx-1,0,-1
        RR(n)=1.d0/((2.d0*n+1.d0)/y-RR(n+1)) ! R(n) := Rn
    enddo
    psim(0)=RR(0)*cos(y)
    do n=1,nstop
        psim(n)=RR(n)*psim(n-1) ! PSI(n) := PSIn
    enddo
    psim(0:nstop)=reshape(psim,(/nstop+1/)) !trim redundant elements
    !---------------------------------------------------------------------------------------------
    

    !-----------Evaluations of partial wave coefficients a and b defined by BH83, Eqs.4.56-4.57---------------------------------------
    do n=1,nstop
      a(n)=((mper_p/(mper_m*m_r)*DD(n)+dble(n)/x)*psi(n)-psi(n-1))/((mper_p/(mper_m*m_r)*DD(n)+dble(n)/x)*qsi(n)-qsi(n-1)) ! BH83, Eq.4.88
      b(n)=((mper_m*m_r/mper_p*DD(n)+dble(n)/x)*psi(n)-psi(n-1))/((mper_m*m_r/mper_p*DD(n)+dble(n)/x)*qsi(n)-qsi(n-1)) ! BH83, Eq.4.88
      c(n)=i*m_r/(psim(n)*(qsi(n-1)-qsi(n)*(mper_m*m_r/mper_p*DD(n)+dble(n)/x)))
      d(n)=i*mper_p/mper_m/(psim(n)*(qsi(n-1)-qsi(n)*(mper_p/(mper_m*m_r)*DD(n)+dble(n)/x)))
    enddo
    !----------------------------------------------------------------------------------------------------------------------------------

    do n=1,nstop
      fn1(n)=(2.d0*n+1.d0)
      fn2(n)=(2.d0*n+1.d0)/(n*(n+1.d0))
    enddo

  !------- Following Qsca,Qext,Qabs,S1,S2 are defined only for nonabsorbing media (Imag(m_m)=0)-------
    Qsca=(2.d0/x**2)*sum(fn1*(abs(a)**2+abs(b)**2)) ! BH83, Eq.4.61
    Qext=(2.d0/x**2)*sum(fn1*real(a+b)) ! BH83, Eq.4.62
    Qabs=Qext-Qsca

    dtheta=pi/dble(nang-1)
    do j=1,nang
      theta=(j-1)*dtheta
      mu=cos(theta)
      pie(1)=1.d0
      pie(2)=3.d0*mu*pie(1)
      tau(1)=mu*pie(1)
      tau(2)=2.d0*mu*pie(2)-3.d0*pie(1)
      do n=3,nstop
          pie(n)=((2.d0*n-1.d0)/(n-1.d0))*mu*pie(n-1)-(dble(n)/(n-1.d0))*pie(n-2)
          tau(n)=n*mu*pie(n)-(n+1.d0)*pie(n-1)
      enddo
      S1(j)=sum(fn2*(a*pie+b*tau)) ! BH83, Eq.4.74
      S2(j)=sum(fn2*(a*tau+b*pie))
    enddo
    !---------------------------------------------------------------------------------------------

    !-------Extinction efficiency based on the optical theorem Qext_opt is defined also for absorbing media--------
    ! extinction efficiency calculated from the optical theorem
    Qext_opt=4.d0/(r_p**2*real(k))*imag(S1(1)*(i/real(k))); ! MISH07 Eq.87
    !--------------------------------------------------------------------------------------------------------------

    !--- Qsca,Qabs,Qext are not defined for absorbing media---
    if(abs(imag(m_m))>1.d-10)then
      Qsca=-1000
      Qabs=-1000
      Qext=-1000
      S1(:)=-1000
      S2(:)=-1000
    endif
    !---------------------------------------------------------

    !write(*,*)'a='
    !do n=1,nstop
      !write(*,"(i0,2e12.4)")n,a(n)
    !enddo
    !write(*,*)'b='
    !do n=1,nstop
      !write(*,"(i0,2e12.4)")n,b(n)
    !enddo
    !write(*,*)'c='
    !do n=1,nstop
      !write(*,"(i0,2e12.4)")n,c(n)
    !enddo
    !write(*,*)'d='
    !do n=1,nstop
      !write(*,"(i0,2e12.4)")n,d(n)
    !enddo

  end subroutine gMIE_core

end module gMIE_module




!==== Following is an example of main routine calling the gMie_core ====

program main
  use gMIE_module
  implicit none

  double precision, parameter :: pi=acos(-1.0d0)

  !---INPUT ARGUMENTS---
  double precision wl_0 ! wavelength in vacuum (=c/w)
  double precision  r_p ! particle radius
  complex(kind(0d0)) eper_m ! complex permittivity ε of medium
  complex(kind(0d0)) mper_m ! complex permiability μ of medium
  complex(kind(0d0)) eper_p ! complex permittivity ε of particle
  complex(kind(0d0)) mper_p ! complex permiability μ of particle
  integer nang ! number of scattering angle between 0-180 deg for scattering matrix element S1 S2
  !---------------------

  double precision Qsca !  scattering efficiency := Csca/(pi*r_p**2)
  double precision Qext ! extinction efficiency := Cext/(pi*r_p**2)
  double precision Qabs ! absorption efficiency := Cabs/(pi*r_p**2)
  double precision Qext_opt ! extinction efficiency applicable to absorbing media := Cext/(pi*r_p**2)
  complex(kind(0d0)),allocatable :: S1(:),S2(:)! element of the 2*2 amplitude scattering matrix defined as BH83, Eq.3.12

  nang=3 !
  allocate(S1(1:nang),S2(1:nang))

  wl_0=1.d0 ! free space wavelength
  r_p=1.d0 ! particle radius
  eper_m=(1.0d0,0.d0) ! complex permittivity ε of medium
  mper_m=1.d0 ! complex permiability μ of medium
  eper_p=(1.d0,0.d0) ! complex permittivity ε of particle
  mper_p=(1.d0,0.d0) ! complex permiability μ of particle

  call gMIE_core(wl_0,r_p,eper_m,mper_m,eper_p,mper_p,nang,Qsca,Qabs,Qext,Qext_opt,S1,S2)

  write(*,'(a6,f10.5)')'Qsca=',Qsca
  write(*,'(a6,f10.5)')'Qabs=',Qabs
  write(*,'(a6,f10.5)')'Qext=',Qext
  write(*,'(a9,f10.5)')'Qext_opt=',Qext_opt

end

!==================================================================
