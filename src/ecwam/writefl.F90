! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,      &
     &                   FILENAME, IUNIT, LOUNIT, LRSTPARAL)

! ---------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!    &              FILENAME,IUNIT,LOUNIT,LRSTPARAL)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *IJINF*     FIRST LOWER DIMENSION OF FL
!     *IJSUP*     FIRST UPPER DIMENSION OF FL
!     *KINF*      SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*      SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*    LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LRSTPARAL* LOGICAL, TRUE THEN WRITE IN PARALLEL

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE (A FORTRAN WRITE)

!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IUNIT.
!     IT IS IN A FORM THAT IS INDEPENDENT OF THE MODEL DECOMPOSITION.



!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE ISO_C_BINDING, ONLY : C_INT, C_LOC, C_PTR

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D, LLUNSTR
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT

      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT), TARGET :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LRSTPARAL

      INTEGER(KIND=JWIM) :: IJ, K, M, J2, J3
      INTEGER(KIND=JWIM) :: KRET, KOUNT, LFILE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), TARGET :: FL_G
      REAL(KIND=JWRB), ALLOCATABLE :: RFL(:,:,:)
      LOGICAL :: LSTREAM

!     Mirror of READFL: swap the whole buffer in one C pass before the
!     existing single WRITE, instead of letting CONVERT='BIG_ENDIAN' swap
!     it element by element. No record markers to emit by hand here either.
      INTEGER(KIND=JWIM), PARAMETER :: IELEMBYTES = STORAGE_SIZE(1.0_JWRB)/8
#if defined(__amdflang__)
      LOGICAL, PARAMETER :: LLFAST = (IELEMBYTES == 4 .OR. IELEMBYTES == 8)
#else
      LOGICAL, PARAMETER :: LLFAST = .FALSE.
#endif
      LOGICAL, PARAMETER :: LLWIDE = (IELEMBYTES == 8)

      INTERFACE
        SUBROUTINE BSWAP32_ARRAY(BUF, N) BIND(C, NAME='bswap32_array')
          USE ISO_C_BINDING, ONLY : C_PTR, C_INT
          TYPE(C_PTR), VALUE :: BUF
          INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: N
        END SUBROUTINE BSWAP32_ARRAY
        SUBROUTINE BSWAP64_ARRAY(BUF, N) BIND(C, NAME='bswap64_array')
          USE ISO_C_BINDING, ONLY : C_PTR, C_INT
          TYPE(C_PTR), VALUE :: BUF
          INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: N
        END SUBROUTINE BSWAP64_ARRAY
      END INTERFACE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WRITEFL',0,ZHOOK_HANDLE)

      LSTREAM = .NOT. LRSTPARAL

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      IF (LLFAST) THEN
!       Bypass IWAM_GET_UNIT's CONVERT='BIG_ENDIAN': the swap happens by
!       hand below, in one pass, instead of element by element on write.
        IF (LOUNIT) THEN
          OPEN(NEWUNIT=IUNIT, FILE=FILENAME(1:LFILE),                   &
     &         FORM='UNFORMATTED', ACCESS='STREAM',                     &
     &         STATUS='REPLACE', ACTION='WRITE')
        ELSE
          OPEN(NEWUNIT=IUNIT, FILE=FILENAME(1:LFILE),                   &
     &         FORM='UNFORMATTED', ACCESS='STREAM', POSITION='APPEND',  &
     &         STATUS='OLD', ACTION='WRITE')
        ENDIF
      ELSEIF (LOUNIT) THEN
        IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0, 'READWRITE', LSTREAM=LSTREAM)
      ELSE
        IUNIT=IWAM_GET_UNIT(-1, FILENAME(1:LFILE) , 'a', 'u', 0, 'READWRITE', LSTREAM=LSTREAM)
      ENDIF

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        FL_G(0,:,:)=0.0_JWRB
        CALL UNBLKRORD(1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,               &
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),             &
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))

        IF (LLFAST) CALL SWAP_ELEMS(C_LOC(FL_G(IJINF,KINF,MINF)), SIZE(FL_G))
        WRITE(IUNIT) (((FL_G(IJ,J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        IF (LLFAST) THEN
!         Swap a copy so the caller's FL is left untouched.
          FL_G = FL
          CALL SWAP_ELEMS(C_LOC(FL_G(IJINF,KINF,MINF)), SIZE(FL_G))
          WRITE(IUNIT) (((FL_G(IJ,J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
        ELSE
          WRITE(IUNIT) (((FL(IJ,J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
        ENDIF
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!     BUT THE SINGLE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(J3, J2, IJ)
        DO J3 = MINF, MSUP
          DO J2 = KINF, KSUP
            DO IJ = IJINF, IJSUP
              FL_G(IJ,J2,J3) = FL(IJ2NEWIJ(IJ),J2,J3)
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        IF (LLFAST) CALL SWAP_ELEMS(C_LOC(FL_G(IJINF,KINF,MINF)), SIZE(FL_G))
        WRITE(IUNIT) FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      ENDIF

      CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('WRITEFL',1,ZHOOK_HANDLE)

      CONTAINS

      SUBROUTINE SWAP_ELEMS(P, N)
!     Byte-reverse N reals at P, at whatever width the build uses. Ported
!     from get_save_spec_opt's byte-swap fast path, with the record-marker
!     handling stripped out: this file format has none to preserve.
      TYPE(C_PTR), INTENT(IN) :: P
      INTEGER(KIND=JWIM), INTENT(IN) :: N
      IF (LLWIDE) THEN
        CALL BSWAP64_ARRAY(P, INT(N, C_INT))
      ELSE
        CALL BSWAP32_ARRAY(P, INT(N, C_INT))
      ENDIF
      END SUBROUTINE SWAP_ELEMS

      END SUBROUTINE WRITEFL
