module io

  public :: RDINPUT
  private

contains

  subroutine RDINPUT(ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, &
       &     PRTEIGVC2, PRTPULS, PRTDIPL, LHWHM, FSTOP, SHIFTP, INFL,  &
       &     TITLE, DIM, FILEAUX, POTCH, POTCHF, POTFILE, POTFILEF,  &
       &     PRTCRL, TDI, TPABSOR, TPCLC, &
       &     TPDIAG, TPPROPG, TPTRANS, TPWIND, TPDIPL, EFC, INIEIGVC, &
       &     DMFILE, IIL, IIU, NIS, IFL, IFU, MF, MP, NFS, NSEED, NSHOT, &
       &     NPR, MXCST, MXDM, NISPG, NFSPG, NREORT, MAXINT, KP, KL, KLX, &
       &     NP, ABSTOL, DT, DTF, TI, TF, TOL, WP, E0, TP, TD, T0, SNI, &
       &     OMG, CHWHM, DE, E1I, ADE, GAMMA, OMEGA, DMX, TPX, T0X, AA, &
       &     QC, DELQ, CSTI, CSTF, CSTFB, CSTDM, SM, SMF, XI, XF, AL, AR, &
       &     rK, rA, X0, VOI, VAR,SHF,PRTREGION,NREG,RANGE,PRTONLYREIM, &
       &     GAMMA3,TDIPOL,OMG3,E03,TP3,TD3,T03,SNI3,KL3)
    IMPLICIT NONE

    LOGICAL  ::        ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, PRTEIGVC2
    LOGICAL  ::        PRTPULS, PRTDIPL, LHWHM, FSTOP, SHIFTP
    LOGICAL  ::        PRTREGION,PRTONLYREIM
    CHARACTER*(*)  ::  INFL, TITLE, DIM, FILEAUX, POTCH, POTCHF, POTFILE
    CHARACTER*(*)  ::  POTFILEF, PRTCRL, TDI, TPABSOR, TPCLC, TPDIAG 
    CHARACTER*(*)  ::  TPPROPG, TPTRANS, TPWIND, TPDIPL, EFC, INIEIGVC 
    CHARACTER*(*)  ::  DMFILE
    INTEGER        ::  IIL, IIU, NIS, IFL, IFU, MF, MP, MXCST, MXDM 
    INTEGER        ::  NISPG, NFS, NFSPG, NSEED, NSHOT, NPR, NREORT
    INTEGER        ::  MAXINT, KP, KL, KLX, NREG
    REAL*8         ::  ABSTOL, DT, DTF, TI, TF, TOL, WP, E0, TP, TD, T0
    REAL*8         ::  SNI, OMG, CHWHM, DE, E1I, ADE
    REAL*8         ::  GAMMA, OMEGA, DMX, TPX, T0X, AA, QC, DELQ
    INTEGER        ::  KL3(3)
    REAL*8         ::  GAMMA3(3),TDIPOL(3)
    real*8                     :: VMIN,VMINB,VMINC
    REAL*8         ::  E03(3),TP3(3),OMG3(3),T03(3),SNI3(3),TD3(3)
    real*8                     :: skl3(3)
    INTEGER        ::  NP(*)
    REAL*8         ::  CSTI(*), CSTF(*), CSTFB(*), SM(*), SMF(*), CSTDM(*)
    REAL*8         ::  XI(*), XF(*), AL(*), AR(*), rK(*), rA(*), X0(*) 
    REAL*8         ::  VOI(*), VAR(*), SHF(*) 
    REAL*8         ::  RANGE(5,7)
    REAL*8,  PARAMETER ::            ZERO     = +0.0D+0, &
         ONE      = +1.0D+0,   &  
         A0A      = +5.291772083D-1, & !A/a.u.
         FATEEVAU = +27.2113834D+0, & !eV/a.u. 
         SHBAR    = +6.58211889D-1 !ev*fs

    INTEGER       I, J,INFLGTH, ITEST, ND
    REAL*8        SKL, SKLX
    REAL*8        WST(MXDM), VFY(10)
    INTEGER       ICHLENGTH, ICOMPAR, IOTEST     
    integer :: ios = 0, line = 0, funit=10, pos ! Note that funit=10 is arbitrary. It could be any unit you want.
    integer :: ierr
    real    :: timescan
    logical :: dim_initialized = .FALSE., seed_initialized=.FALSE., nshot_initialized=.FALSE.
    character(len=128) :: buffer, another_buffer, label    
    character(len=32)  :: filename


    call get_command_argument(1, filename)
    print *, " * Reading the input file : ", filename
    open(funit, file=filename, iostat=ierr); call chkioerr(ierr,info="Reading input file")
    ! Iterate over the keywords in the filename for the first time and get non-array parameters
    do while (ios == 0) ! while file still has lines to read
       read(funit, '(A)', iostat=ios) buffer ! we read the line and put it in the buffer
       if (ios == 0) then  ! And if ios == 0, there is a line. If it's not, file has ended.
          select case (buffer)     ! Now let's iterate over all possible labels and read them
          case ('*TITLE')      
             read(funit, *, iostat=ierr) TITLE; call chkioerr(ierr,info="*TITLE")
             write(*,*) " Title = ", TITLE
          case ('*TPCALC')      
             read(funit, *, iostat=ierr) TPCLC; call chkioerr(ierr,info="*TPCALC")
             write(*,*) " TPCALC = ", TPCLC
             select case (TPCLC)
             case ('.ENERGY')
                TDI='.TI'
             case ('.PROPAGATION','.CORRELATION','.COLLISION')
                TDI='.TD'
             case('.SPECTRUM')
                read(funit, *, iostat=ierr) TDI; call chkioerr(ierr,info=".SPECTRUM: TDI")
             case('.ONLYSPEC')
                read(funit, *, iostat=ierr) FILEAUX; call chkioerr(ierr,info=".ONLYSPEC: FILEAUX")
                TDI = '.NC'
             case('.PES','.POTENTIAL')
                TDI='.NC'
             case default
                write(*,*) "TPCLC = ", TPCLC, "; Exiting, do not know how to proceed. "
                stop
             end select
          case ('*DIMENSION')
             read(funit, *, iostat=ierr) DIM; call chkioerr(ierr,info="*DIMENSION")
             write(*,*) " [*DIMENSION]: ", DIM
             select case (DIM)
             case ('.1D')
                read(funit, *, iostat=ierr) NP(1); call chkioerr(ierr,info="*DIMENSION: NP(1)")
             case ('.2D')
                read(funit, *, iostat=ierr) NP(1),NP(2); call chkioerr(ierr,info="*DIMENSION: NP(1),NP(2)")
             case('.3D')
                read(funit, *, iostat=ierr) NP(1),NP(2),NP(3); call chkioerr(ierr,info="*DIMENSION: NP(1),NP(2),NP(3)")
             case default
                write(*,*) "DIM = ", DIM, "; Exiting, do not know how to proceed. "
                stop
             end select
             dim_initialized=.TRUE.
          case ('*GRID_CUT')
             read(funit, *, iostat=ierr) (VAR(I),I=1,MXCST,1); call chkioerr(ierr,info="*GRID_CUT")
          case('*TPTRANS')
             read(funit, *, iostat=ierr) TPTRANS; call chkioerr(ierr,info="*TPTRANS")
          case('*NREORT')
             read(funit, *, iostat=ierr) NREORT; call chkioerr(ierr,info="*NREORT")
          case('*FINAL_POTENTIAL')
             read(funit, *, iostat=ierr) POTCH; call chkioerr(ierr,info="*FINAL_POTENTIAL")
             select case(POTCH)
             case ('.FILE')
                read(funit, *, iostat=ierr) POTFILE; call chkioerr(ierr,info="*FINAL_POTENTIAL: POTFILE")
             case default
                read(funit, *, iostat=ierr) (CSTF(I),I=1,MXCST,1); call chkioerr(ierr,info="*FINAL_POTENTIAL: CSTF")
             end select
          case('*POTENTIAL')
             read(funit, *, iostat=ierr) POTCH; call chkioerr(ierr,info="*POTENTIAL")
             write(*,*) " Reading potential from: ", POTCH
             POTCHF = POTCH
             select case(POTCHF)
             case('.FILE')
                read(funit, *, iostat=ierr) POTFILE; call chkioerr(ierr,info="*POTENTIAL: POTFILEF")
                POTFILEF=POTFILE
             case default
                read(funit, *, iostat=ierr) (CSTI(I),I=1,MXCST,1); call chkioerr(ierr,info="*POTENTIAL: CSTI")
                read(funit, *, iostat=ierr) (CSTF(I),I=1,MXCST,1); call chkioerr(ierr,info="*POTENTIAL: CSTF")
             end select
          case('*GRID_RANGES')
             read(funit, *, iostat=ierr) (XI(I), XF(I), I=1,MXDM,1); call chkioerr(ierr,info="*GRID_RANGES")
             write(*,*) " *GRID_RANGES: ", xi(1), xf(1), xi(2), xf(2)
             if (buffer(14:15).eq.'au') then ! WTF?
                DO I=1,MXDM,1
                   XI(I) = XI(I)*A0A
                   XF(I) = XF(I)*A0A
                   write(*,*)xi(i),xf(i)
                ENDDO
             end if
          case('*INIEIGVC')
             read(funit, *, iostat=ierr) INIEIGVC; call chkioerr(ierr,info="*INIEIGVC")
          case('*MASS_FACTOR')
             read(funit, *, iostat=ierr) (SMF(I),I=1,MXCST,1); call chkioerr(ierr,info="*MASS_FACTOR")
          case('*MASS')
             read(funit, *, iostat=ierr) (SM(I),I=1,MXCST,1); call chkioerr(ierr,info="*MASS")
          case('*SCALING_FACTOR')
             read(funit, *, iostat=ierr) (SHF(I),I=1,MXCST,1); call chkioerr(ierr,info="*SCALING_FACTOR")
          case('*CROSS_TERM')
             read(funit, *, iostat=ierr) SHF(3); call chkioerr(ierr,info="*CROSS_TERM") ! WTF? SHF was filled at *SCALING_FACTOR.
          case('*SEED')
             if (seed_initialized) then
                write(*,*) " ERROR: *SEED was already initialized. Check the input. "
                stop
             end if
             read(funit, *, iostat=ierr) NSEED; call chkioerr(ierr,info="*SEED"); seed_initialized=.TRUE.
          case('*NSHOT')
             if (nshot_initialized) then
                write(*,*) " ERROR: *NSHOT was already initialized. Check the input. "
                stop
             end if
             read(funit, *, iostat=ierr) NSHOT; call chkioerr(ierr,info="*NSHOT"); nshot_initialized=.TRUE. 
          case('*CHANGE')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*CHANGE")
             select case (another_buffer)
             case('.NO','.OFF')
                CHANGE = .FALSE.
             case default  ! Bad practice.
                CHANGE = .TRUE.
             end select
          case('*SHIFT_POTENTIAL')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*SHIFT_POTENTIAL")
             select case (another_buffer)
             case('.NO','.OFF')
                SHIFTP = .FALSE.
             case default ! Bad practice again.
                SHIFTP = .TRUE.
             end select
          case('*PRTEIGVC2')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTEIGVC2")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTEIGVC2= .TRUE.
             case default
                PRTEIGVC2= .FALSE.
             end select
          case('*PRTREGION')
             read(funit, *, iostat=ierr) NREG ; call chkioerr(ierr,info="*PRTREGION")
             if (dim_initialized) then
                select case (DIM)
                case ('.1D')
                   do i=1,Nreg
                      read(funit, *, iostat=ierr) (RANGE(I,J), J=1,2,1) ; call chkioerr(ierr,info="*PRTREGION: 1D")
                   end do

                case ('.2D')
                   do i=1,Nreg
                      read(funit, *, iostat=ierr) (RANGE(I,J), J=1,4,1) ; call chkioerr(ierr,info="*PRTREGION: 2D")
                   end do
                case ('.3D')
                   do i=1,Nreg
                      read(funit, *, iostat=ierr) (RANGE(I,J), J=1,6,1) ; call chkioerr(ierr,info="*PRTREGION: 3D")
                   end do
                end select
             else
                write(*,*) " ERROR: DIM must be initialized before *PRTREGION. Exiting."
                stop
             end if

          case ('*PRTONLYREIM')
             PRTONLYREIM = .TRUE.
          case ('*PRTEIGVC')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTEIGVC")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTEGVC = .TRUE.
             case default
                PRTEGVC = .FALSE.
             end select
          case ('*PRTPOT')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTPOT")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTPOT = .TRUE.
             case default
                PRTPOT = .FALSE.
             end select
          case ('*PRTVEFF')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTVEFF")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTVEFF = .TRUE.
             case default
                PRTVEFF = .FALSE.
             end select
          case ('*PRTPULSE')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTPULSE")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTPULS = .TRUE.
             case default
                PRTPULS = .FALSE.
             end select
          case ('*PRTDIPOLE')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTDIPOLE")
             select case (another_buffer)
             case ('.YES','.ON')
                PRTDIPL = .TRUE.
             case default
                PRTDIPL = .FALSE.
             end select
          case ('*DONOTSTOP')
             FSTOP = read_logical(funit, info='*DONOTSTOP')
          case ('*PRTCRL')
             read(funit, *, iostat=ierr) another_buffer ; call chkioerr(ierr,info="*PRTCRL")
             select case (another_buffer)
             case ('.FULL','.ALL')
                PRTCRL='.YES'
             case ('.NONE','.NULL')
                PRTCRL='.NO'
             end select
          case ('*INITIAL_POTENTIAL','*INITIAL_CONDITION','*INITIAL_WP','*INITIAL_WF',&
               '*INITIAL_WAVE_FUNCTION','*WAVE_FUNCTION','*WAVE_PACKET')
             write(*,*) " ERROR: Using one of the wrong keywrods."
             write(*,*) " To feed in the potenital, use the *POTENTIAL keyword. Potentials must be located in one file. "
             stop
             ! read(funit, *, iostat=ierr) POTCH ; call chkioerr(ierr,info="*(Initial potential, 8 keywords)")
             ! select case (POTCH)
             ! case ('.FILE')
             !    read(funit, *, iostat=ierr) POTFILE ; call chkioerr(ierr,info="*(Initial potential, 8 keywords) : .File")
             !    write(*,*) " [*POTENTIAL] Potential filename = ", POTFILE
             ! case default
             !    read(funit, *, iostat=ierr) (CSTI(I), I=1,MXCST,1) ; call chkioerr(ierr,info="*(Initial potential, 8 keywords) : CSTI(:)")
             ! end select
          case ('*TPDIAG')
             read(funit, *, iostat=ierr) TPDIAG ; call chkioerr(ierr,info="*TPDIAG")
          case ('*NIST')
             read(funit, *, iostat=ierr) NIS, IIL ; call chkioerr(ierr,info="*NIST")
             IIU = IIL + NIS - 1
          case ('*NFST')
             read(funit, *, iostat=ierr) NFS, IFL ; call chkioerr(ierr,info="*NFST")
             IFU = IFL + NFS - 1
          case ('*ABSTOL')
             read(funit, *, iostat=ierr) ABSTOL ; call chkioerr(ierr,info="*ABSTOL")
          case ('*MAXINT')
             read(funit, *, iostat=ierr) MAXINT ; call chkioerr(ierr,info="*MAXINT")
          case ('*COLLISION_ENERGY')
             read(funit, *, iostat=ierr) (rK(i),i=1,MXDM); call chkioerr(ierr,info="*COLLISION_ENERGY")
             IF(buffer(19:20).NE.'au')THEN
                DO I=1,MXDM,1
                   rK(I) = rK(I)/FATEEVAU 
                ENDDO
             ENDIF
          case ('*INIT_POSSITION')
             read(funit, *, iostat=ierr) (X0(i),i=1,MXDM); call chkioerr(ierr,info="*INIT_POSSITION")
          case ('*SIZE_OF_GAUSSIAN')
             read(funit, *, iostat=ierr) (rA(i),i=1,MXDM); call chkioerr(ierr,info="*SIZE_OF_GAUSSIAN")
             IF(buffer(19:20).EQ.'au')THEN
                DO I=1,MXDM,1
                   rA(i)=rA(i)*A0A
                ENDDO
             ENDIF

          case ('*PROPAG')
             read(funit, *, iostat=ierr) TPPROPG; call chkioerr(ierr,info="*PROPAG")
             select case (TPPROPG)
             case ('.PSIL','.PPSIL','.PPLNZ','.PSH','.PPSH')
                read(funit, *, iostat=ierr) TI, TF, DT, MP; call chkioerr(ierr,info="*PROPAG: TPPROPG (.PSIL ...)")
             case ('.PSOD','.PFSOD','.PWSOD','.PPSOD','.PFSPO','.PWSPO') 
                read(funit, *, iostat=ierr) TI, TF, DT; call chkioerr(ierr,info="*PROPAG: (.PSOD ...)")
             case ('.P2P')
                read(funit, *, iostat=ierr) TI, TF, DT; call chkioerr(ierr,info="*PROPAG: (.P2P)")
                read(funit, *, iostat=ierr) another_buffer; call chkioerr(ierr,info="*PROPAG: (.P2P: next line)")
                select case (another_buffer)
                case ('.PULSE')
                   read(funit, *, iostat=ierr) TPX, T0X, SKLX; call chkioerr(ierr,info="*PROPAG: (.P2P: .PULSE)")
                   KLX = SKLX
                   VFY(1) = ZERO
                case ('.GAMMA_OMEG_TRNSDIPOL')
                   read(funit, *, iostat=ierr) GAMMA, OMEGA, DMX; call chkioerr(ierr,info="*PROPAG: (.P2P: .GAMMA_OMEG_TRNSDIPOL)")
                   GAMMA = GAMMA/FATEEVAU
                   OMEGA = OMEGA/SHBAR
                   DMX   = DMX/FATEEVAU
                   VFY(2) = ZERO
                case ('.AA_QC_DELQ')
                   read(funit, *, iostat=ierr) AA, QC, DELQ; call chkioerr(ierr,info="*PROPAG: (.P2P: .AA_QC_DELQ)")
                   if (buffer(12:15).ne.('(au)')) then
                      AA = AA/FATEEVAU
                      QC = QC/A0A
                      DELQ = DELQ/A0A
                   end if
                   VFY(3) = ZERO
                end select
             case ('.COUPLE3')
                read(funit, *, iostat=ierr) TI, TF, DT; call chkioerr(ierr,info="*PROPAG: (.COUPLE3)")
                write(*,*) " [.COUPLE3]: TI, TF, DT = ", TI, TF, DT
                !read(funit, *, iostat=ierr) another_buffer; call chkioerr(ierr,info="*PROPAG: (.COUPLE3: next line)")
             end select
          case('.DECAYS')
             read(funit, *, iostat=ierr) (GAMMA3(i), i=1,3); call chkioerr(ierr,info="*PROPAG: (.COUPLE3: .DECAYS)")
             write(*,*) " GAMMA3 = ", GAMMA3(:)
          case('.PULSES')
             read(funit, *, iostat=ierr) E03(1), TP3(1), OMG3(1), T03(1), SNI3(1), SKL3(1);
             call chkioerr(ierr,info="*PROPAG: (.COUPLE3: .PULSES (1) )")
             read(funit, *, iostat=ierr) E03(2), TP3(2), OMG3(2), T03(2), SNI3(2), SKL3(2);
             call chkioerr(ierr,info="*PROPAG: (.COUPLE3: .PULSES (2) )")
             KL3(1) = SKL3(1)
             KL3(2) = SKL3(2)
             write(*,*) " * Pulses output: " 
             write(*,*) " E03 = ", E03(:)
             write(*,*) " TP3 = ", TP3(:)
             write(*,*) " OMG3 = ", OMG3(:)
             write(*,*) " T03 = ", T03(:)
             write(*,*) " SNI3 = ", SNI3(:)
             write(*,*) " SKL3 = ", SKL3(:)
             write(*,*) " * "
          case ('.TDIPOL')
             read(funit, *, iostat=ierr) TDIPOL(1), TDIPOL(2); call chkioerr(ierr,info="*PROPAG: (.COUPLE3: .TDIPOL)")
          case ('*ABSORBING')
             read(funit, *, iostat=ierr) TPABSOR; call chkioerr(ierr,info="*ABSORBING")
             select case (TPABSOR)
             case ('.SMOOTHW','.SMRS','.VOPTIC')
                ABSORB = .TRUE.
                read(funit, *, iostat=ierr) (VOI(I), I=1,MXDM,1); call chkioerr(ierr,info="*ABSORBING: (.SMOOTHW: VOI)")
                read(funit, *, iostat=ierr) (AL(I), AR(I), I=1,MXDM,1); call chkioerr(ierr,info="*ABSORBING: (.SMOOTHW: AL, AR)")
             case ('.NULL')
                ABSORB = .FALSE.
                ! WTF? Reading in the same variable.
                read(funit, *, iostat=ierr) (WST(I), I=1,MXDM,1); call chkioerr(ierr,info="*ABSORBING: (.NULL: WST)")
                read(funit, *, iostat=ierr) (WST(I), WST(I), I=1,MXDM,1); call chkioerr(ierr,info="*ABSORBING: (.NULL: WST WST)")
             case default
                write(*,*) " ERROR: wrong keywords for *ABSORBING. Check the input. "
                stop
             end select
          case ('*PULSE_AND_DIPOLE')
             read(funit, *, iostat=ierr) EFC; call chkioerr(ierr,info="*PULSE_AND_DIPOLE")
             select case (EFC)
             case ('.GAUS')
                read(funit, *, iostat=ierr) E0, TP, OMG, T0, SNI; call chkioerr(ierr,info="*PULSE_AND_DIPOLE: .GAUS")
             case ('.GGAUS')
                read(funit, *, iostat=ierr) E0, TP, OMG, T0, SNI, SKL; call chkioerr(ierr,info="*PULSE_AND_DIPOLE: .GGAUS")
             case ('.SIN2')
                read(funit, *, iostat=ierr) E0, TP, OMG, TD, SNI; call chkioerr(ierr,info="*PULSE_AND_DIPOLE: .SIN2")
             case ('.NONE')
                read(funit, *, iostat=ierr) (WST(I),I=1,MXDM); call chkioerr(ierr,info="*PULSE_AND_DIPOLE: .NONE")
             case default
                write(*,*) " ERROR: *PULSE_AND_DIPOLE. Unknown keyword. Check the input"
                stop
             end select
             read(funit, *, iostat=ierr) TPDIPL; call chkioerr(ierr,info="*PULSE_AND_DIPOLE: TPDIPL")
             select case (TPDIPL)
             case ('.READ','.GET')
                read(funit, *, iostat=ierr) DMFILE; call chkioerr(ierr,info="*PULSE_AND_DIPOLE: (TPDIPL: .READ) ")
             case ('.NULL')
                CSTDM(1) = ZERO
             case default
                read(funit, *, iostat=ierr) (CSTDM(i), i=1,MXCST); call chkioerr(ierr,info="*PULSE_AND_DIPOLE: (TPDIPL: default) ")
             end select
          case ('*PRPGSTATE')
             read(funit, *, iostat=ierr) NISPG; call chkioerr(ierr,info="*PRPGSTATE ")
          case ('*PRPTOL')
             read(funit, *, iostat=ierr) TOL; call chkioerr(ierr,info="*PRPTOL ")
          case ('*NPROJECTIONS')
             read(funit, *, iostat=ierr) NPR, KP; call chkioerr(ierr,info="*NPROJECTIONS ")
          case ('*FOURIER')
             read(funit, *, iostat=ierr) MF; call chkioerr(ierr,info="*FOURIER")
          case ('*WINDOWING')
             read(funit, *, iostat=ierr) TPWIND; call chkioerr(ierr,info="*WINDOWING: TPWIND")
             read(funit, *, iostat=ierr) WP; call chkioerr(ierr,info="*WINDOWING: WP")
          case ('*HWHM')
             read(funit, *, iostat=ierr) CHWHM, DE, E1I, ADE; call chkioerr(ierr,info="*HWHM: CHWHM...")
             LHWHM = .TRUE.

          end select

       end if
    end do

    ! ### Arbitrary code used in the RDINPUT.f ### 
    IF(DIM(1:3).EQ.'.1D')THEN
       DO I=2,MXDM,1
          NP(I) = ONE
       ENDDO
    ELSEIF(DIM(1:3).EQ.'.2D')THEN
       DO I=3,MXDM,1
          NP(I) = ONE
       ENDDO
    ELSEIF(TPCLC(1:9).EQ.'.ONLYSPEC')THEN
       CONTINUE
    ELSE
       WRITE(*,*) "Unknown error: DIM = "
       WRITE(*,*)DIM
       IF(FSTOP) STOP
    ENDIF

    IF(INIEIGVC(1:4).EQ.'.GET' .OR. INIEIGVC(1:5).EQ.'.READ'&
       &     .OR. TPCLC(1:9).EQ.'.ONLYSPEC')THEN
       NIS = 1
    ELSEIF(INIEIGVC(1:5).EQ.'.CALC')THEN
       CONTINUE
    ELSE
       WRITE(*,*) "Unknown error: INIEIGVC = "
       WRITE(*,*)INIEIGVC
       IF(FSTOP) STOP
    ENDIF

    IF(NPR+KP.GT.NIS)NPR = NPR - NIS
    IF(TPPROPG(1:7).EQ.'.P2PABM' .AND. VFY(1) .NE. ZERO .AND.&
       &     TPPROPG(1:8).NE.'.P2PABM2')THEN
       WRITE(*,*)'It is necessary specify the pulse parameters!'
       IF(FSTOP) STOP
    ELSEIF(TPPROPG(1:7).EQ.'.P2PABM' .AND. VFY(2) .NE. ZERO .AND.&
       &        TPPROPG(1:8).NE.'.P2PABM2')THEN
       WRITE(*,*)'<<<>>> It is necessary specify the decay rate, frequ&
       &ency and transition the dipole moment! <<<>>>'
       IF(FSTOP) STOP
    ENDIF

    IF(TPPROPG(1:7).EQ.'.P2PSOD' .AND. VFY(3) .NE. ZERO .OR.&
       &     TPPROPG(1:8).EQ.'.P2PABM2' .AND. VFY(3) .NE. ZERO)THEN
       WRITE(*,*)'<<<>>> It is necessary specify the coupling constant&
       &s! <<<>>>'
       IF(FSTOP) STOP
    ENDIF

  end subroutine RDINPUT

  logical function read_logical(funit, info)
    integer, intent(in) :: funit
    character(128) :: buffer
    character(*), intent(in), optional   :: info
    integer :: ierr
    read(funit, *, iostat=ierr) buffer 
    if (present(info)) then
       call chkioerr(ierr,info=info)
    else
       call chkioerr(ierr)
    end if
    
    select case (buffer)
    case ('.YES','.ON')
       read_logical = .TRUE.
    case default
       read_logical = .FALSE.
    end select
  end function read_logical
  
  
  subroutine chkioerr(ierr, info)
    integer, intent(in) :: ierr
    character(len=*), intent(in), optional :: info
    if (ierr.ne.0) then
       write(*,"(A, I4)") " ERROR: I/O failed with ierr=", ierr
       if (present(info)) then
          write(*,*) " Info: ", info
       end if
       stop
    end if
  end subroutine chkioerr

end module io
