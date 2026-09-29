!> LFP half-cell model: Li foil | separator | porous LiFePO4 cathode.
!>
!> A self-contained program: edit the input file, run, and it writes Time_Voltage.txt.
!> The equations are described in docs/model.md; parameters in docs/parameters.md.
!> The linear block system of each time step is solved with bandsolver (Newman's BAND).
!>
!>   usage:  lfp_f [input.nml]        (default input file: lfp.nml)
!>
!> mode = 'faithful' reproduces the original research code, including the defects
!> listed in docs/deviations.md; mode = 'corrected' applies the fixes.
!>
!> SPDX-License-Identifier: BSD-3-Clause
program lfp
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use bandsolver_kernel, only: band_solve, BAND_OK, PIVOT_LEGACY
    implicit none

    integer, parameter :: NV = 4                  ! unknowns per node
    integer, parameter :: IC = 1, IP1 = 2, IP2 = 3, ICS = 4

    ! ---------------- parameters (defaults = original research code) ----------------
    real(dp) :: L_cath_um = 24.0_dp, L_sep = 25.0e-4_dp     ! cathode thickness in um (L_cath = L_cath_um*1e-4 cm)
    integer  :: nj = 101, sep_node = 22
    real(dp) :: eps = 0.5_dp, eps_AM = 0.8_dp, eps_sep = 0.39_dp, tau_sep = 4.0_dp, bruggeman = -0.5_dp
    real(dp) :: f_AM = 0.8_dp          ! corrected mode: active fraction of the solid phase (D-9)
    real(dp) :: D = 2.0e-6_dp, t_plus = 0.25_dp, c_bulk = 1.0e-3_dp, z_plus = 1.0_dp, z_minus = -1.0_dp
    real(dp) :: sigma = 3.0e-3_dp, M = 125.759_dp, rho = 3.6_dp, Q_th = 0.170_dp, R_p = 200.0e-7_dp
    real(dp) :: k_rxn = -1.0_dp                    ! < 0: use the default for the chosen mode
    real(dp) :: alpha_a = 0.5_dp, alpha_c = 0.5_dp, k_Li = 1.0e-6_dp, c_Li_ref = 1.0e-3_dp
    real(dp) :: R = 8.314_dp, T = 298.0_dp, F = 96485.0_dp
    real(dp) :: C_rate = 1.0_dp, phi1_init = 3.6_dp, phi2_init = 0.0_dp, cs_init = 1.0e-5_dp
    real(dp) :: t_max = 36000.0_dp
    integer  :: n_steps = 36000
    real(dp) :: V_min = 2.5_dp, V_max = 4.2_dp        ! corrected-mode cutoffs (D-3)
    real(dp) :: fd_step = 1.0e-6_dp
    real(dp) :: newton_tol = 1.0e-10_dp                ! corrected-mode Newton tolerance (D-7)
    integer  :: newton_max_iter = 25
    character(len=16)  :: mode = 'faithful'
    character(len=2048) :: steps = ''                  ! corrected-mode protocol (docs/protocol.md)
    integer  :: cycles = 1
    real(dp) :: write_interval = 18.0_dp               ! [s], corrected mode
    character(len=256) :: file = 'Time_Voltage.txt'

    namelist /cell/ L_cath_um, L_sep, nj, sep_node, eps, eps_AM, f_AM, eps_sep, tau_sep, bruggeman
    namelist /electrolyte/ D, t_plus, c_bulk, z_plus, z_minus
    namelist /active/ sigma, M, rho, Q_th, R_p, k_rxn, alpha_a, alpha_c, k_Li, c_Li_ref
    namelist /constants/ R, T, F
    namelist /operation/ C_rate, phi1_init, phi2_init, cs_init, t_max, n_steps, V_min, V_max
    namelist /numerics/ fd_step, newton_tol, newton_max_iter, mode
    namelist /protocol/ steps, cycles
    namelist /output/ file, write_interval

    ! ---------------- derived quantities and mode switches ----------------
    logical  :: faithful
    real(dp) :: lit36, ocp_c(7), spec_a, tortuosity, i_spec, i_app, eps_sep_face, phi1_sign
    logical  :: full_current
    real(dp) :: dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c
    real(dp) :: L_cath                   ! cathode thickness [cm], from L_cath_um as in the original (24 * 1.0d-4)
    real(dp) :: mass_area, i_1C, vf_AM   ! vf_AM: active volume fraction (eps_AM or f_AM*(1-eps))
    real(dp), parameter :: THETA_REG = 1.0e-6_dp       ! D-13 regularization threshold

    ! ---------------- protocol (corrected mode) ----------------
    integer, parameter :: K_CC = 1, K_CV = 2, K_REST = 3
    integer :: nstep
    integer, allocatable :: skind(:)
    real(dp), allocatable :: sC(:), sV(:), sT(:), sVmin(:), sVmax(:), sImin(:)   ! sT, sImin < 0: unset

    ! ---------------- mesh and state ----------------
    integer :: s                                  ! interface node index (1-based = sep_node)
    real(dp), allocatable :: dx(:), aW(:), aE(:), bW(:), bE(:)
    real(dp), allocatable :: c(:,:), dc(:,:), A(:,:,:), B(:,:,:), Dm(:,:,:), G(:,:)

    ! ---------------- time loop ----------------
    integer  :: ounit, nsolve
    real(dp) :: time, dt, mAhg
    character(len=1) :: state
    character(len=64) :: exit_reason
    character(len=256) :: input_file

    ! ---------------- read input ----------------
    input_file = 'lfp.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    call read_input(trim(input_file))
    call setup()

    open(newunit=ounit, file=trim(file), status='replace', action='write')

    allocate(c(NV,nj), dc(NV,nj), A(NV,NV,nj), B(NV,NV,nj), Dm(NV,NV,nj), G(NV,nj))
    c(IC,:) = c_bulk
    c(IP1,:) = phi1_init
    c(IP2,:) = phi2_init
    c(ICS,:) = cs_init
    dc = 0.0_dp
    nsolve = 0
    time = 0.0_dp
    mAhg = 0.0_dp
    dt = t_max/real(n_steps, dp)
    exit_reason = 'max_steps'

    if (faithful) then
        call run_faithful()
    else
        call parse_protocol()
        call run_protocol()
    end if

    close(ounit)
    write(*,'(A,G0,A,I0,A)') trim(mode)//' run, C-rate ', C_rate, ': exit '//trim(exit_reason)//' after ', &
        nsolve, ' steps; wrote '//trim(file)

contains

    ! =============================== faithful run ===============================
    subroutine run_faithful()
        !! The original program's constant-current discharge, step for step.
        integer :: it, status
        integer :: last_write                 ! an integer in the original (D-4)
        real(dp) :: write_every
        last_write = 0
        write_every = t_max/n_steps/200
        state = 'D'
        if (C_rate < 0) state = 'C'
        do it = 1, n_steps
            if (it == 1) then
                call write_row(.true.)
            else if ((time - last_write)/3600 >= write_every) then
                call write_row(.false.)
                last_write = int(time - dt)
            else if (it >= n_steps) then
                call write_row(.false.)
            else if (c(IP1,nj) >= 99.0_dp .and. state == 'C') then
                call write_row(.false.)
                exit_reason = 'end_of_charge'
                exit
            else if (ieee_is_nan(dc(IC,1))) then
                call write_row(.false.)
                exit_reason = 'nan'
                exit
            else if (time >= 99.0_dp*3600.0_dp) then
                call write_row(.false.)
                exit_reason = 'max_time'
                exit
            end if
            if (state == 'D') then
                mAhg = mAhg + 1000.0_dp*i_spec*dt/3600.0_dp
            else if (state == 'C') then
                mAhg = mAhg - 1000.0_dp*i_spec*dt/3600.0_dp
            end if
            call assemble(dt)
            call band_solve(NV, nj, A, B, Dm, G, dc, status, pivot=PIVOT_LEGACY)
            if (status /= BAND_OK) dc = ieee_value(1.0_dp, ieee_quiet_nan)
            c = c + dc
            nsolve = it
            if (state == 'R') then
                dt = dt*1.0001_dp
            else
                dt = t_max/real(n_steps, dp)
            end if
            time = time + dt
        end do
    end subroutine run_faithful

    ! =============================== protocol (corrected mode) ===============================
    subroutine parse_protocol()
        !! Parse `steps` (docs/protocol.md) and expand it `cycles` times.
        character(len=len(steps)) :: txt, part
        character(len=64) :: word, key
        integer :: n1, k, pos, semi, cyc, w0, w1, eqp, ios
        real(dp) :: val
        real(dp), allocatable :: tC(:), tV(:), tT(:), tVmin(:), tVmax(:), tImin(:)
        integer, allocatable :: tk(:)
        logical :: hasC, hasV

        txt = adjustl(steps)
        if (len_trim(txt) == 0) then
            allocate(skind(1), sC(1), sV(1), sT(1), sVmin(1), sVmax(1), sImin(1))
            nstep = 1
            skind = K_CC; sC = C_rate; sV = 0.0_dp; sT = -1.0_dp; sVmin = V_min; sVmax = V_max; sImin = -1.0_dp
            return
        end if
        n1 = count_steps(txt)
        allocate(tk(n1), tC(n1), tV(n1), tT(n1), tVmin(n1), tVmax(n1), tImin(n1))
        pos = 1
        k = 0
        do while (pos <= len_trim(txt))
            semi = index(txt(pos:), ';')
            if (semi == 0) then
                part = txt(pos:)
                pos = len_trim(txt) + 1
            else
                part = txt(pos:pos+semi-2)
                pos = pos + semi
            end if
            part = adjustl(part)
            if (len_trim(part) == 0) cycle
            k = k + 1
            tC(k) = 0.0_dp; tV(k) = 0.0_dp; tT(k) = -1.0_dp; tVmin(k) = V_min; tVmax(k) = V_max; tImin(k) = -1.0_dp
            hasC = .false.; hasV = .false.
            w0 = 1
            call next_word(part, w0, w1, word)
            select case (lower(word))
            case ('cc');   tk(k) = K_CC
            case ('cv');   tk(k) = K_CV
            case ('rest'); tk(k) = K_REST
            case default;  call proto_error(k, 'unknown step type '//trim(word))
            end select
            w0 = w1
            do
                call next_word(part, w0, w1, word)
                if (len_trim(word) == 0) exit
                w0 = w1
                eqp = index(word, '=')
                if (eqp == 0) call proto_error(k, 'expected key=value, got '//trim(word))
                key = lower(word(:eqp-1))
                word = word(eqp+1:)
                call fix_exponent(word)
                read(word, *, iostat=ios) val
                if (ios /= 0) call proto_error(k, 'bad number '//trim(word))
                select case (tk(k))
                case (K_CC)
                    select case (key)
                    case ('c');    tC(k) = val; hasC = .true.
                    case ('t');    tT(k) = val
                    case ('vmin'); tVmin(k) = val
                    case ('vmax'); tVmax(k) = val
                    case default;  call proto_error(k, 'cc does not take '//trim(key))
                    end select
                case (K_CV)
                    select case (key)
                    case ('v');    tV(k) = val; hasV = .true.
                    case ('t');    tT(k) = val
                    case ('imin'); tImin(k) = val
                    case default;  call proto_error(k, 'cv does not take '//trim(key))
                    end select
                case (K_REST)
                    if (key /= 't') call proto_error(k, 'rest does not take '//trim(key))
                    tT(k) = val
                end select
                if (key == 't' .and. val <= 0) call proto_error(k, 't must be positive')
            end do
            if (tk(k) == K_CC .and. .not. hasC) call proto_error(k, 'cc needs C=')
            if (tk(k) == K_CV .and. .not. hasV) call proto_error(k, 'cv needs V=')
            if (tk(k) == K_CV .and. tT(k) < 0 .and. tImin(k) < 0) call proto_error(k, 'cv needs t= or Imin= to end')
            if (tk(k) == K_REST .and. tT(k) < 0) call proto_error(k, 'rest needs t=')
        end do
        cyc = max(1, cycles)
        nstep = k*cyc
        allocate(skind(nstep), sC(nstep), sV(nstep), sT(nstep), sVmin(nstep), sVmax(nstep), sImin(nstep))
        skind = [(tk(1:k), n1 = 1, cyc)]
        sC = [(tC(1:k), n1 = 1, cyc)]
        sV = [(tV(1:k), n1 = 1, cyc)]
        sT = [(tT(1:k), n1 = 1, cyc)]
        sVmin = [(tVmin(1:k), n1 = 1, cyc)]
        sVmax = [(tVmax(1:k), n1 = 1, cyc)]
        sImin = [(tImin(1:k), n1 = 1, cyc)]
    end subroutine parse_protocol

    integer function count_steps(txt)
        character(len=*), intent(in) :: txt
        integer :: i
        count_steps = 1
        do i = 1, len_trim(txt)
            if (txt(i:i) == ';') count_steps = count_steps + 1
        end do
    end function count_steps

    subroutine next_word(str, i0, i1, word)
        !! The next blank-separated word of str at or after position i0; i1 is the position after it.
        character(len=*), intent(in) :: str
        integer, intent(in) :: i0
        integer, intent(out) :: i1
        character(len=*), intent(out) :: word
        integer :: i, j
        word = ''
        i = i0
        do while (i <= len_trim(str))
            if (str(i:i) /= ' ') exit
            i = i + 1
        end do
        j = i
        do while (j <= len_trim(str))
            if (str(j:j) == ' ') exit
            j = j + 1
        end do
        if (j > i) word = str(i:j-1)
        i1 = j
    end subroutine next_word

    function lower(str) result(out)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: out
        integer :: i
        out = str
        do i = 1, len(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') out(i:i) = achar(iachar(str(i:i)) + 32)
        end do
    end function lower

    subroutine fix_exponent(word)
        !! Accept Fortran d exponents written in lower case as well.
        character(len=*), intent(inout) :: word
        integer :: i
        do i = 1, len_trim(word)
            if (word(i:i) == 'd' .or. word(i:i) == 'D') word(i:i) = 'e'
        end do
    end subroutine fix_exponent

    subroutine proto_error(k, msg)
        integer, intent(in) :: k
        character(len=*), intent(in) :: msg
        write(error_unit,'(A,I0,A)') 'protocol step ', k, ': '//msg
        error stop 2
    end subroutine proto_error

    subroutine run_protocol()
        !! Run the protocol steps in order (docs/protocol.md).
        integer :: k, nsteps_done
        real(dp) :: I, h, h_done, t_step, last_write
        logical :: stopped, ok
        character(len=16) :: why, reason

        I = 0.0_dp
        if (skind(1) == K_CC) I = sC(1)*i_1C
        i_app = I
        call write_row_c(.true., 1)
        last_write = time
        nsteps_done = 0
        reason = ''
        do k = 1, nstep
            t_step = 0.0_dp
            if (skind(k) == K_CC) then
                I = sC(k)*i_1C
            else if (skind(k) == K_REST) then
                I = 0.0_dp
            end if
            do
                h = dt
                if (sT(k) >= 0) h = min(dt, sT(k) - t_step)
                if (skind(k) == K_CV) then
                    call cv_step(h, sV(k), I, ok)
                    h_done = h
                    stopped = sImin(k) >= 0 .and. abs(I) <= sImin(k)*i_1C
                    why = 'current_limit'
                else
                    i_app = I
                    call advance(h, h_done, stopped, ok, skind(k) == K_CC, sVmin(k), sVmax(k))
                    why = 'cutoff_high'
                    if (stopped) then
                        if (cell_voltage() <= sVmin(k)) why = 'cutoff_low'
                    end if
                end if
                i_app = I
                if (.not. ok) then
                    call write_row_c(.false., k)
                    exit_reason = 'solver_fail'
                    nsolve = nsteps_done
                    return
                end if
                mAhg = mAhg + 1000.0_dp*(I/mass_area)*h_done/3600.0_dp
                time = time + h_done
                t_step = t_step + h_done
                nsteps_done = nsteps_done + 1
                if (any(ieee_is_nan(c))) then
                    call write_row_c(.false., k)
                    exit_reason = 'nan'
                    nsolve = nsteps_done
                    return
                end if
                if (stopped .or. (sT(k) >= 0 .and. t_step >= sT(k)*(1.0_dp - 1.0e-12_dp))) then
                    call write_row_c(.false., k)
                    last_write = time
                    if (stopped) then
                        reason = why
                    else
                        reason = 'duration'
                    end if
                    exit
                end if
                if (time - last_write >= write_interval) then
                    call write_row_c(.false., k)
                    last_write = time
                end if
                if (time >= 99.0_dp*3600.0_dp) then
                    call write_row_c(.false., k)
                    exit_reason = 'max_time'
                    nsolve = nsteps_done
                    return
                end if
            end do
        end do
        if (nstep == 1) then
            exit_reason = reason
        else
            exit_reason = 'end_of_protocol'
        end if
        nsolve = nsteps_done
    end subroutine run_protocol

    subroutine cv_step(h, V_set, I, ok)
        !! One constant-voltage time step: find I with V(I) = V_set (see simulate.cv_step in Python).
        real(dp), intent(in) :: h, V_set
        real(dp), intent(inout) :: I
        logical, intent(out) :: ok
        real(dp), parameter :: tol = 1.0e-9_dp
        real(dp) :: c_start(NV,nj), fI, grow, a, b, fa, fb
        logical :: have_a, have_b, good
        integer :: it, side
        c_start = c
        ok = .false.
        call cv_feval(h, V_set, c_start, I, fI, good)
        if (good .and. abs(fI) <= tol) then
            ok = .true.
            return
        end if
        grow = max(abs(I), 1.0e-2_dp*i_1C)
        have_a = .false.; have_b = .false.
        a = 0; b = 0; fa = 0; fb = 0
        do it = 1, 60
            if (fI > 0) then
                a = I; fa = fI; have_a = .true.
                if (have_b) exit
                if (I < 0) then
                    I = 0.0_dp
                else
                    I = I + grow
                end if
            else
                b = I; fb = fI; have_b = .true.
                if (have_a) exit
                if (I > 0) then
                    I = 0.0_dp
                else
                    I = I - grow
                end if
            end if
            grow = grow*2.0_dp
            call cv_feval(h, V_set, c_start, I, fI, good)
            if (good .and. abs(fI) <= tol) then
                ok = .true.
                return
            end if
        end do
        if (.not. (have_a .and. have_b)) then
            c = c_start
            return
        end if
        side = 0
        do it = 1, 200
            if (ieee_is_finite(fa) .and. ieee_is_finite(fb)) then
                I = (a*fb - b*fa)/(fb - fa)
                if (.not. (a < I .and. I < b)) I = 0.5_dp*(a + b)
            else
                I = 0.5_dp*(a + b)
            end if
            call cv_feval(h, V_set, c_start, I, fI, good)
            if (abs(fI) <= tol .or. (b - a) <= 1.0e-14_dp*i_1C) then
                ok = good
                if (.not. good) c = c_start
                return
            end if
            if (fI > 0) then
                a = I; fa = fI
                if (side == 1 .and. ieee_is_finite(fb)) fb = fb*0.5_dp
                side = 1
            else
                b = I; fb = fI
                if (side == -1 .and. ieee_is_finite(fa)) fa = fa*0.5_dp
                side = -1
            end if
        end do
        c = c_start
    end subroutine cv_step

    subroutine cv_feval(h, V_set, c_start, Itry, fval, success)
        !! f(I) = V - V_set after one Newton step from c_start; +/-inf if the step fails.
        real(dp), intent(in) :: h, V_set, c_start(NV,nj), Itry
        real(dp), intent(out) :: fval
        logical, intent(out) :: success
        c = c_start
        i_app = Itry
        call newton_step(h, success)
        if (success) then
            fval = cell_voltage() - V_set
        else
            fval = ieee_value(1.0_dp, ieee_positive_inf)
            if (Itry >= 0) fval = -fval
        end if
    end subroutine cv_feval

    ! =============================== input ===============================
    subroutine read_input(path)
        character(len=*), intent(in) :: path
        integer :: u, ios
        logical :: exists
        inquire(file=path, exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') 'input file not found: '//path
            error stop 2
        end if
        open(newunit=u, file=path, status='old', action='read')
        rewind(u); read(u, nml=cell, iostat=ios);        call check(ios, 'cell')
        rewind(u); read(u, nml=electrolyte, iostat=ios); call check(ios, 'electrolyte')
        rewind(u); read(u, nml=active, iostat=ios);      call check(ios, 'active')
        rewind(u); read(u, nml=constants, iostat=ios);   call check(ios, 'constants')
        rewind(u); read(u, nml=operation, iostat=ios);   call check(ios, 'operation')
        rewind(u); read(u, nml=numerics, iostat=ios);    call check(ios, 'numerics')
        rewind(u); read(u, nml=protocol, iostat=ios);    call check(ios, 'protocol')
        rewind(u); read(u, nml=output, iostat=ios);      call check(ios, 'output')
        close(u)
    end subroutine read_input

    subroutine check(ios, group)
        integer, intent(in) :: ios
        character(len=*), intent(in) :: group
        if (ios > 0) then
            write(error_unit,'(A)') 'error reading namelist group &'//group
            error stop 2
        end if
        ! ios < 0: group absent, keep defaults
    end subroutine check

    pure real(dp) function r32(x)
        !! Round to single precision (reproduces the original's single-precision literals).
        real(dp), intent(in) :: x
        r32 = real(real(x, sp), dp)
    end function r32

    ! =============================== setup ===============================
    subroutine setup()
        real(dp), parameter :: OCP_FIT(7) = [3.114559_dp, 4.438792_dp, 71.7352_dp, 70.85337_dp, &
                                             4.240252_dp, 68.5605_dp, 67.730082_dp]
        ! 10**0.966 as the original evaluated it: single precision, constant-folded
        real(dp), parameter :: K_FAITHFUL = 1.0e-8_dp*real(10.0_sp**0.966_sp, dp)
        real(dp) :: h_sep, h_cat, t_an, d_cat, d_an, u_cat, u_an
        integer :: j

        select case (trim(mode))
        case ('faithful')
            faithful = .true.
        case ('corrected')
            faithful = .false.
        case default
            write(error_unit,'(A)') 'mode must be ''faithful'' or ''corrected'''
            error stop 2
        end select

        if (faithful) then
            R = r32(R); c_bulk = r32(c_bulk); Q_th = r32(Q_th); M = r32(M); rho = r32(rho)
            phi1_init = r32(phi1_init); eps_sep = r32(eps_sep); eps_AM = r32(eps_AM); C_rate = r32(C_rate)
            if (k_rxn < 0) k_rxn = K_FAITHFUL
            lit36 = r32(3.6_dp)
            ocp_c = [(r32(OCP_FIT(j)), j = 1, 7)]
            eps_sep_face = eps          ! D-2
            phi1_sign = -1.0_dp         ! D-1
            full_current = .false.      ! D-11
        else
            if (k_rxn < 0) k_rxn = 1.0e-8_dp*10.0_dp**0.966_dp
            lit36 = 3.6_dp
            ocp_c = OCP_FIT
            eps_sep_face = eps_sep
            phi1_sign = 1.0_dp
            full_current = .true.
        end if

        if (faithful) then
            vf_AM = eps_AM
        else
            vf_AM = f_AM*(1.0_dp - eps)
        end if
        L_cath = L_cath_um*1.0e-4_dp
        spec_a = 3*vf_AM/R_p
        tortuosity = eps**bruggeman
        i_spec = Q_th*C_rate
        i_app = i_spec*L_cath*vf_AM*rho
        mass_area = L_cath*vf_AM*rho
        i_1C = Q_th*mass_area

        ! ion diffusivities and mobilities, then effective values per region
        t_an = 1.0_dp - t_plus
        d_cat = D*(1.0_dp + (t_an/t_plus))/(2.0_dp*t_an/t_plus)
        d_an = d_cat*t_an/t_plus
        u_cat = d_cat/(R*T)
        u_an = d_an/(R*T)
        dcat_s = d_cat/tau_sep;    dan_s = d_an/tau_sep;    ucat_s = u_cat/tau_sep;    uan_s = u_an/tau_sep
        dcat_c = d_cat/tortuosity; dan_c = d_an/tortuosity; ucat_c = u_cat/tortuosity; uan_c = u_an/tortuosity

        s = sep_node
        allocate(dx(nj), aW(nj), aE(nj), bW(nj), bE(nj))
        h_sep = L_sep/real(sep_node - 2, dp)
        h_cat = L_cath/real(nj - sep_node - 1, dp)
        dx = 0.0_dp
        dx(2:s-1) = h_sep
        dx(s+1:nj-1) = h_cat
        aW = 0.0_dp; bW = 0.0_dp; aE = 0.0_dp; bE = 0.0_dp
        do j = 2, nj
            aW(j) = dx(j-1)/(dx(j-1) + dx(j))
            bW(j) = 2.0_dp/(dx(j-1) + dx(j))
        end do
        do j = 1, nj - 1
            aE(j) = dx(j)/(dx(j+1) + dx(j))
            bE(j) = 2.0_dp/(dx(j) + dx(j+1))
        end do
    end subroutine setup

    ! =============================== kinetics ===============================
    real(dp) function cs_max()
        cs_max = (rho/M)*M*Q_th*1000.0_dp*lit36/F
    end function cs_max

    real(dp) function ocp(cs, cc)
        !! Open-circuit potential [V]; corrected mode adds the Nernst term (RT/F) ln(c/c_bulk) (D-14).
        real(dp), intent(in) :: cs, cc
        real(dp) :: th
        th = (cs/(rho/M))/(M*Q_th*1000.0_dp*lit36/F)
        ocp = ocp_c(1) + ocp_c(2)*atan(-(ocp_c(3)*th) + ocp_c(4)) - ocp_c(5)*atan(-(ocp_c(6)*th) + ocp_c(7))
        if (.not. faithful) ocp = ocp + R*T/F*log(cc/c_bulk)
    end function ocp

    real(dp) function rate(cc, cs, p1, p2)
        !! Butler-Volmer current density per interfacial area [A/cm2], anodic positive.
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp) :: eta, i0
        eta = p1 - p2 - ocp(cs, cc)
        i0 = F*k_rxn*(cc**alpha_a)*((cs_max() - cs)**alpha_a)*(cs**alpha_c)
        if (faithful) i0 = r32(i0)     ! D-4: implicitly single precision in the original
        rate = i0*(exp(alpha_a*F*eta/(R*T)) - exp(-(alpha_c*F*eta/(R*T))))
    end function rate

    subroutine rate_derivs(cc, cs, p1, p2, i, di)
        !! Rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: h
        if (.not. faithful) then
            call rate_derivs_exact(cc, cs, p1, p2, i, di)
            return
        end if
        h = fd_step
        i = rate(cc, cs, p1, p2)
        if (cc <= h) then
            di(IC) = (rate(cc + h, cs, p1, p2) - i)/h
        else
            di(IC) = (rate(cc + h, cs, p1, p2) - rate(cc - h, cs, p1, p2))/(2.0_dp*h)
        end if
        if (cs <= h) then
            di(ICS) = (rate(cc, cs + h, p1, p2) - i)/h
        else
            di(ICS) = (rate(cc, cs + h, p1, p2) - rate(cc, cs - h, p1, p2))/(2.0_dp*h)
        end if
        di(IP1) = (rate(cc, cs, p1 + h, p2) - rate(cc, cs, p1 - h, p2))/(2.0_dp*h)
        di(IP2) = (rate(cc, cs, p1, p2 + h) - rate(cc, cs, p1, p2 - h))/(2.0_dp*h)
    end subroutine rate_derivs

    real(dp) function ocp_slope(cs)
        !! dU/dcs [V cm3/mol]
        real(dp), intent(in) :: cs
        real(dp) :: th, dth, du
        th = (cs/(rho/M))/(M*Q_th*1000.0_dp*lit36/F)
        dth = 1.0_dp/((rho/M)*(M*Q_th*1000.0_dp*lit36/F))
        du = ocp_c(2)*(-ocp_c(3))/(1.0_dp + (-(ocp_c(3)*th) + ocp_c(4))**2) &
           - ocp_c(5)*(-ocp_c(6))/(1.0_dp + (-(ocp_c(6)*th) + ocp_c(7))**2)
        ocp_slope = du*dth
    end function ocp_slope

    subroutine power_reg(x, alpha, delta, g, dg)
        !! x**alpha, replaced below delta by a C1 quadratic with g(0) = 0 and a finite slope (D-13).
        real(dp), intent(in) :: x, alpha, delta
        real(dp), intent(out) :: g, dg
        real(dp) :: u
        if (x < delta) then
            u = x/delta
            g = delta**alpha*((2.0_dp - alpha)*u + (alpha - 1.0_dp)*u*u)
            dg = delta**(alpha - 1.0_dp)*((2.0_dp - alpha) + 2.0_dp*(alpha - 1.0_dp)*u)
        else
            g = x**alpha
            dg = alpha*x**(alpha - 1.0_dp)
        end if
    end subroutine power_reg

    subroutine rate_derivs_exact(cc, cs, p1, p2, i, di)
        !! Rate and exact derivatives w.r.t. (c, phi1, phi2, cs); corrected mode (fixes D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: rt, aa, bb, eta, i0, ea, ec, di_deta, gv, dgv, gs, dgs, pre, di0
        rt = R*T
        aa = alpha_a*F/rt
        bb = alpha_c*F/rt
        eta = p1 - p2 - ocp(cs, cc)
        call power_reg(cs_max() - cs, alpha_a, THETA_REG*cs_max(), gv, dgv)
        call power_reg(cs, alpha_c, THETA_REG*cs_max(), gs, dgs)
        pre = F*k_rxn*(cc**alpha_a)
        i0 = pre*gv*gs
        di0 = pre*(gs*(-dgv) + gv*dgs)
        ea = exp(aa*eta)
        ec = exp(-bb*eta)
        i = i0*(ea - ec)
        di_deta = i0*(aa*ea + bb*ec)
        di(IC) = alpha_a*i/cc - di_deta*(rt/F)/cc      ! includes dU/dc of the Nernst term
        di(ICS) = di0*(ea - ec) - di_deta*ocp_slope(cs)
        di(IP1) = di_deta
        di(IP2) = -di_deta
    end subroutine rate_derivs_exact

    ! =============================== corrected-mode time step ===============================
    subroutine time_terms(dt, Tt)
        !! Storage coefficients: the time-derivative part of each row is Tt*(c - c_old).
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: Tt(NV,nj)
        Tt = 0.0_dp
        Tt(ICS,1) = -(vf_AM/dt)
        Tt(IC,2:s-1) = -(eps_sep/dt*dx(2:s-1))
        Tt(ICS,2:s-1) = -((1.0_dp - eps_sep)/dt)
        Tt(IC,s+1:nj-1) = -((eps/dt)*dx(s+1:nj-1))
        Tt(ICS,s:nj) = -(vf_AM/dt)
    end subroutine time_terms

    real(dp) function bounded_step(dcc)
        !! Largest step <= 1 keeping 0 < c and 0 < cs < cs_max (at most 90 % of the way to a bound)
        !! and changing no potential by more than 0.1 V (Butler-Volmer exponentials make Newton overshoot).
        real(dp), intent(in) :: dcc(NV,nj)
        real(dp), parameter :: keep = 0.9_dp
        real(dp) :: csm
        integer :: j
        real(dp), parameter :: max_dphi = 0.1_dp   ! largest potential change per iteration [V]
        real(dp) :: dphi
        csm = cs_max()
        bounded_step = 1.0_dp
        dphi = maxval(abs(dcc(IP1:IP2,:)))
        if (dphi > max_dphi) bounded_step = max_dphi/dphi
        do j = 1, nj
            if (dcc(IC,j) < 0) bounded_step = min(bounded_step, keep*(c(IC,j) - 0.0_dp)/(-dcc(IC,j)))
            if (dcc(ICS,j) < 0) bounded_step = min(bounded_step, keep*(c(ICS,j) - 0.0_dp)/(-dcc(ICS,j)))
            if (dcc(ICS,j) > 0) bounded_step = min(bounded_step, keep*(csm - c(ICS,j))/dcc(ICS,j))
        end do
        bounded_step = max(bounded_step, 0.0_dp)
    end function bounded_step

    subroutine equilibrate()
        !! Scale every equation by the largest entry of its row in B (the solution is unchanged).
        integer :: j, r
        real(dp) :: sc
        do j = 1, nj
            do r = 1, NV
                sc = maxval(abs(B(r,:,j)))
                if (sc == 0.0_dp) sc = 1.0_dp
                A(r,:,j) = A(r,:,j)/sc
                B(r,:,j) = B(r,:,j)/sc
                Dm(r,:,j) = Dm(r,:,j)/sc
                G(r,j) = G(r,j)/sc
            end do
        end do
    end subroutine equilibrate

    subroutine newton_step(h, ok)
        !! One backward-Euler step of length h, solved with Newton's method; c is updated on success.
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        real(dp) :: c_old(NV,nj), Tt(NV,nj), scale(NV), lam
        integer :: k, st
        c_old = c
        call time_terms(h, Tt)
        scale = [c_bulk, 1.0_dp, 1.0_dp, cs_max()]
        ok = .false.
        do k = 1, newton_max_iter
            call assemble(h)
            G = G - Tt*(c - c_old)
            call equilibrate()
            call band_solve(NV, nj, A, B, Dm, G, dc, st)
            if (st /= BAND_OK) exit
            lam = bounded_step(dc)
            c = c + lam*dc
            if (lam == 1.0_dp .and. maxval(abs(dc)/spread(scale, 2, nj)) <= newton_tol) then
                ok = .true.
                return
            end if
        end do
        c = c_old
    end subroutine newton_step

    subroutine advance(dt, t_done, stopped, ok, check, vlo, vhi)
        !! Advance by dt at the current i_app, halving the sub-step on Newton failure. With `check`,
        !! a sub-step that crosses a voltage cutoff by more than 0.1 mV is halved, so the step ends
        !! within 0.1 mV of the cutoff (see simulate.advance in Python).
        real(dp), intent(in) :: dt, vlo, vhi
        logical, intent(in) :: check
        real(dp), intent(out) :: t_done
        logical, intent(out) :: stopped, ok
        real(dp), parameter :: min_dt = 1.0e-6_dp, event_dv = 1.0e-4_dp, event_min_dt = 1.0e-9_dp
        real(dp) :: hh, vv, mg, c_save(NV,nj)
        logical :: good
        t_done = 0.0_dp
        hh = dt
        stopped = .false.
        ok = .true.
        do while (t_done < dt)
            hh = min(hh, dt - t_done)
            c_save = c
            call newton_step(hh, good)
            if (.not. good) then
                if (hh/2 < min_dt) then
                    ok = .false.
                    return
                end if
                hh = hh/2
                cycle
            end if
            if (check) then
                vv = cell_voltage()
                mg = min(vv - vlo, vhi - vv)
                if (mg < 0.0_dp) then
                    if (mg < -event_dv .and. hh/2 >= event_min_dt) then
                        c = c_save
                        hh = hh/2
                        cycle
                    end if
                    t_done = t_done + hh
                    stopped = .true.
                    return
                end if
            end if
            t_done = t_done + hh
        end do
    end subroutine advance

    ! =============================== assembly ===============================
    subroutine face_coeffs(e, dcat, ucat, dan, uan, cface, gphi2, dd, ff)
        !! Cation-flux (row 1) and ionic-current (row 3) coefficients of one face.
        real(dp), intent(in) :: e, dcat, ucat, dan, uan, cface, gphi2
        real(dp), intent(inout) :: dd(NV,NV), ff(NV,NV)
        real(dp) :: k
        k = z_plus**2*ucat + z_minus**2*uan
        dd(IC,IC)   = -(e*dcat)
        ff(IC,IC)   = -(e*z_plus*ucat*F*gphi2)
        dd(IC,IP2)  = -(e*z_plus*ucat*F*cface)
        dd(IP2,IC)  = -(e*F*(z_plus*dcat + z_minus*dan))
        ff(IP2,IC)  = -(e*F**2*k*gphi2)
        dd(IP2,IP2) = -(e*F**2*k*cface)
    end subroutine face_coeffs

    subroutine assemble(dt)
        real(dp), intent(in) :: dt
        real(dp) :: dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), rj(NV,NV), gg(NV)
        real(dp) :: cW(NV), cE(NV), gW(NV), gE(NV), i, di(NV)
        integer :: j, k

        A = 0.0_dp; B = 0.0_dp; Dm = 0.0_dp; G = 0.0_dp
        do j = 1, nj
            dW = 0.0_dp; dE = 0.0_dp; fW = 0.0_dp; fE = 0.0_dp; rj = 0.0_dp; gg = 0.0_dp
            cW = 0.0_dp; cE = 0.0_dp; gW = 0.0_dp; gE = 0.0_dp
            if (j > 1) then
                cW = aW(j)*c(:,j) + (1.0_dp - aW(j))*c(:,j-1)
                gW = bW(j)*(c(:,j) - c(:,j-1))
            end if
            if (j < nj) then
                cE = aE(j)*c(:,j+1) + (1.0_dp - aE(j))*c(:,j)
                gE = bE(j)*(c(:,j+1) - c(:,j))
            end if
            call rate_derivs(c(IC,j), c(ICS,j), c(IP1,j), c(IP2,j), i, di)

            if (j == 1) then
                ! ---- Li-foil face ----
                dE(IC,IC)  = -(eps_sep_face*dcat_s)
                fE(IC,IC)  = -(eps_sep_face*z_plus*ucat_s*F*gE(IP2))
                dE(IC,IP2) = -(eps_sep_face*z_plus*ucat_s*F*cE(IC))
                gg(IC) = -i_app/F + (dE(IC,IC)*gE(IC) + fE(IC,IC)*cE(IC))
                rj(ICS,ICS) = 0.0_dp - 1.0_dp*vf_AM/dt
                dE(IP1,IP1) = -(1.0_dp - eps_sep_face)*sigma
                gg(IP1) = phi1_sign*dE(IP1,IP1)*gE(IP1)
                rj(IP2,IP2) = 1.0_dp
                gg(IP2) = 0.0_dp - c(IP2,j)

            else if (j < s) then
                ! ---- separator interior ----
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cE(IC), gE(IP2), dE, fE)
                rj(IC,IC) = -(eps_sep/dt*dx(j))
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                rj(ICS,ICS) = -(1.0_dp*(1.0_dp - eps_sep)/dt)
                dW(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)

            else if (j == s) then
                ! ---- separator/cathode interface ----
                call face_coeffs(eps_sep_face, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps_sep_face)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else if (j < nj) then
                ! ---- cathode interior ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                do k = 1, NV
                    rj(IC,k) = (spec_a*di(k)/F)*dx(j)
                    rj(IP1,k) = -((spec_a*di(k))*dx(j))
                    rj(IP2,k) = (spec_a*di(k))*dx(j)
                end do
                rj(IC,IC) = (spec_a*di(IC)/F)*dx(j) - (eps/dt)*dx(j)
                gg(IC) = -((spec_a*i/F)*dx(j)) - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) &
                                               + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = (spec_a*i)*dx(j) - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                           + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = -((spec_a*i)*dx(j)) - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else
                ! ---- current collector ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                gg(IC) = 0.0_dp - dW(IC,IC)*gW(IC) - fW(IC,IC)*cW(IC)
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = i_app - dW(IP1,IP1)*gW(IP1)
                gg(IP2) = 0.0_dp - dW(IP2,IC)*gW(IC) - fW(IP2,IC)*cW(IC)
                call solid_row(i, di, dt, rj, gg)
            end if

            ! control-volume blocks: A dc(j-1) + B dc(j) + D dc(j+1) = G
            if (j == 1) then
                B(:,:,j) = rj - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            else if (j == nj) then
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW
            else
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            end if
            G(:,j) = gg
        end do
    end subroutine assemble

    real(dp) function current(ff, dd, cf, gf)
        !! Ionic-current face value used in the residual (D-11: faithful mode omits diffusion).
        real(dp), intent(in) :: ff(NV,NV), dd(NV,NV), cf(NV), gf(NV)
        current = ff(IP2,IP2)*cf(IP2) + dd(IP2,IP2)*gf(IP2)
        if (full_current) current = current + dd(IP2,IC)*gf(IC)
    end function current

    subroutine solid_row(i, di, dt, rj, gg)
        !! vf_AM dcs/dt = -a i_n / F  (uniform particles)
        real(dp), intent(in) :: i, di(NV), dt
        real(dp), intent(inout) :: rj(NV,NV), gg(NV)
        integer :: k
        do k = 1, NV
            rj(ICS,k) = -(spec_a*di(k)/F)
        end do
        rj(ICS,ICS) = -(spec_a*di(ICS)/F) - 1.0_dp*vf_AM/dt
        gg(ICS) = +(spec_a*i/F)
    end subroutine solid_row

    ! =============================== output ===============================
    real(dp) function li_eta()
        !! Signed overpotential of the lithium counter electrode (negative on discharge).
        real(dp) :: i0_li, alpha
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        alpha = 0.5_dp
        if (.not. faithful) then           ! symmetric Butler-Volmer (D-12); i_app is the present current
            li_eta = -(R*T/(alpha*F))*asinh(i_app/(2.0_dp*i0_li))
        else if (state == 'C') then
            li_eta = 0.5_dp*log(i_app/i0_li)/(alpha*F/(R*T))
        else if (state == 'D') then
            li_eta = -(0.5_dp*log(i_app/i0_li))/(alpha*F/(R*T))
        else
            li_eta = 0.0_dp
        end if
    end function li_eta

    real(dp) function li_nernst()
        !! Corrected mode: Nernst potential of the lithium foil, (RT/F) ln(c/c_Li_ref) [V].
        li_nernst = R*T/F*log(c(IC,1)/c_Li_ref)
    end function li_nernst

    real(dp) function cell_voltage()
        !! Cell voltage against the lithium foil (0 V). The solver fixes the gauge with phi2 = 0 at the
        !! foil face; the equations depend only on potential differences, so the foil-referenced
        !! potentials are the solved ones minus U_Li + eta_Li (the foil metal on the solver's scale).
        !! (Imposing the foil reference as the boundary condition makes the node-1 block singular at rest.)
        if (faithful) then
            cell_voltage = c(IP1,nj) + li_eta()
        else
            cell_voltage = c(IP1,nj) + li_eta() - li_nernst()
        end if
    end function cell_voltage

    subroutine write_row_c(header, step)
        !! Corrected-mode output row: the original columns plus the current and the step index.
        logical, intent(in) :: header
        integer, intent(in) :: step
        real(dp) :: i0_li, eta
        character(len=1) :: st
        if (header) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0', 'Current', 'Step', 'Li_Nernst'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3', 'mA/cm2', '#', 'mV'
        end if
        st = 'R'
        if (i_app > 0) st = 'D'
        if (i_app < 0) st = 'C'
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        eta = li_eta()
        write(ounit,'(A5,1X,2(F12.5,1X),5(ES15.5,1X),I15,1X,ES15.5)') st, time/3600.0_dp, cell_voltage(), &
            mAhg*M*3.6_dp/F, eta*1.0e3_dp, i0_li*1.0e3_dp, c(IC,1), i_app*1.0e3_dp, step, li_nernst()*1.0e3_dp
    end subroutine write_row_c

    subroutine write_row(header)
        logical, intent(in) :: header
        real(dp) :: equiv, i0_li, eta
        if (header) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3'
        end if
        equiv = mAhg*M*lit36/F
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        eta = li_eta()
        write(ounit,'(A5,1X,2(F12.5,1X),20(ES15.5,1X))') state, time/real(3600, dp), c(IP1,nj) + eta, equiv, &
            eta*1.0e3_dp, i0_li*1.0e3_dp, c(IC,1)
    end subroutine write_row

end program lfp
