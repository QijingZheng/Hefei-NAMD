! precision used in vasp

module prec
  integer, parameter :: q =SELECTED_real_KIND(10)
  integer, parameter :: qs=SELECTED_real_KIND(5)
end module


module constants
  use prec

  complex(q), parameter :: imgUnit = (0.0_q, 1.0_q)
  complex(q), parameter :: cero = (0.0_q, 0.0_q)
  complex(q), parameter :: uno = (1.0_q, 0.0_q)
  complex(q), parameter :: hbar = 0.6582119281559802_q ! reduced Planck constant eV / fs
  !  Some important Parameters, to convert to a.u.  !  - AUTOA  = 1. a.u. in Angstroem
  !  - RYTOEV = 1 Ry in Ev
  !  - EVTOJ  = 1 eV in Joule
  !  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
  !  - BOLKEV = Boltzmanns constant in eV/K
  !  - BOLK   = Boltzmanns constant in Joule/K
  real(q), parameter :: AUTOA=0.529177249_q, RYTOEV=13.605826_q
  real(q), parameter :: CLIGHT = 137.037  ! speed of light in a.u.
  real(q), parameter :: EVTOJ=1.60217733E-19_q, AMTOKG=1.6605402E-27_q, &
                        BOLKEV=8.6173857E-5_q, BOLK=BOLKEV*EVTOJ
  real(q), parameter :: EVTOKCAL=23.06

  ! FELECT = (the electronic charge)/(4*pi*the permittivity of free space)
  !         in atomic units this is just e^2
  ! EDEPS = electron charge divided by the permittivity of free space
  !         in atomic units this is just 4 pi e^2
  ! HSQDTM = (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
  real(q), parameter  :: PI =3.141592653589793238_q, TPI=2*PI
  real(q), parameter  :: FELECT = 2*AUTOA*RYTOEV, EDEPS=4*PI*2*RYTOEV*AUTOA, &
                         HSQDTM = RYTOEV*AUTOA*AUTOA
  complex(q), parameter  :: CITPI = (0._q,1._q)*TPI

  ! vector field A times momentum times e/ (2 m_e c) is an energy
  ! magnetic moments are supplied in Bohr magnetons
  ! e / (2 m_e c) A(r) p(r) = energy
  ! e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla =
  ! e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3 = energy
  ! conversion factor from magnetic moment to energy
  ! checked independently in SI by Gilles de Wijs
  real(q), parameter :: MAGMOMTOENERGY=1/CLIGHT**2*AUTOA**3*RYTOEV

  ! dimensionless number connecting input and output magnetic moments
  ! AUTOA e^2 (2 m_e c^2)
  real(q), parameter :: MOMTOMOM=AUTOA/CLIGHT/CLIGHT/2
  real(q), parameter :: AUTOA2=AUTOA *AUTOA
  real(q), parameter :: AUTOA3=AUTOA2*AUTOA
  real(q), parameter :: AUTOA4=AUTOA2*AUTOA2
  real(q), parameter :: AUTOA5=AUTOA3*AUTOA2
  
  ! dipole moment in atomic units to Debye
  real(q), parameter :: AUTDEBYE=2.541746

end module constants
