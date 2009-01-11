C *************************************************************
      SUBROUTINE PSREPHB (EPOCH, ITELNO, RA50, DEC50, PBEPOCH, PB, 
     &                       PBDOT, EPBIN, PBIN, ASINI, WBIN, ECC, 
     &                       TOBS, POBS, XMA, BTDB)
C *************************************************************

      IMPLICIT NONE

      REAL*8 EPOCH, RA50, DEC50, PBEPOCH, PB, PBDOT, EPBIN, PBIN, 
     &       ASINI, WBIN, ECC, TOBS, POBS, XMA, BTDB, RA20, DEC20

      INTEGER ITELNO

*+
* Subroutine to compute the arrival time and period at the telescope.
* All arguments are in double precision.
* Taken from ARTHUR on 13 Mar 1990. The only change is the source of the 
*   telescope information.
* Renamed from PSREPH with extra dummy argument for barycentric TDB.
* Also entry PSREPHBJ for J2000 coordinates. RNM  3 April, 1991.

*
* Given:
*    EPOCH       (DP sec)  The MJD of the observation.
*    ITELNO      (I)       Telescope number.
*    RA50,DEC50  (DP deg)  The B1950 RA and DEC of the pulsar in degrees.
*    RA20,DEC20  (DP deg)  The J2000 RA and DEC of the pulsar in degrees.
*    PBEPOCH     (DP sec)  The barycentric MJD of a pulse and epoch of PB.
*    PB          (DP sec)  The barycentric period at PBEPOCH. If negative,
*                            PBEPOCH, PB and PBDOT are already topocentric.
*    PBDOT       (DP s/s)  The barycentric period derivative.
*    EPBIN       (DP sec)  The barycentric MJD of binary periastron.
*    PBIN        (DP sec)  The binary barycentric period.
*    ASINI       (DP sec)  The binary major semi-axis.
*    WBIN        (DP deg)  The binary longitude of periastron.
*    ECC         (DP  - )  The binary orbital eccentricity.
*
* Returned:
*    TOBS        (DP sec)  The MJD of the next pulse at the telescope.
*    POBS        (DP sec)  The period of the pulses at the telescope.
*    XMA         (DP deg)  Mean anomaly.
*    BTDB        (DP sec)  Equivalent barycentric TDB MJD of observation
*-

* Include telescope information.
      INCLUDE 'tel.def'

* Define astronomical and mathematical constants.
      REAL*8 AUS, ETTAI, DMCONST, TWOPI, RTOD, STRAD, SARAD, REARTH,
     &       C, SOLSID
      PARAMETER ( AUS     = 499.0047837D0
     &           ,ETTAI   = 32.184D0
     &           ,DMCONST = 2.4104616D-4
     &           ,TWOPI   = 6.2831853072D0
     &           ,RTOD  = 360D0/TWOPI
     &           ,STRAD   = TWOPI/86400D0
     &           ,SARAD   = STRAD/15D0
     &           ,REARTH  = 6378.16D0
     &           ,C       = 299792.5D0
     &           ,SOLSID  = 1.00273790931D0 )

      REAL*8 sla_DTT, sla_GMST, sla_DVDV, ARRTIME

* Local variables.
      REAL*8 R1950, D1950, R2000, D2000, DVB(3), DPB(3), DVH(3), DPH(3),
     &       DPS(3), DPVE(6), OBSLONG, OBSLAT, OBSHT, EVEL, EDELAY,
     &       EACCN, BVEL, TDT, ST, BDELAY, BACCN, PTDB, PTIME, 
     &       TDB, DOPP

      REAL*4 G

      INTEGER I

      if (itelno.gt.ntel) then
         write (*,*) "psrephb (Barycentric corrections) error:",
     +        "unknown telescope. Exiting."
         stop
      end if
* Convert coordinates to radians.
      R1950 = RA50/RTOD
      D1950 = DEC50/RTOD

* Convert source coords to J2000.0
      CALL sla_FK45Z (R1950, D1950, 1950.0D0, R2000, D2000)

      GO TO 10

* Entry for J2000 coordinates
      ENTRY PSREPHBJ (EPOCH, ITELNO, RA20, DEC20, PBEPOCH, PB, 
     &                       PBDOT, EPBIN, PBIN, ASINI, WBIN, ECC, 
     &                       TOBS, POBS, XMA, BTDB)
      if (itelno.gt.ntel) then
         write (*,*) "psrephbj (Barycentric corrections) error:",
     +        "unknown telescope. Exiting."
         stop
      end if

* Convert to radians
      R2000 = RA20/RTOD
      D2000 = DEC20/RTOD

* Convert the time to TDT.
10    TDT = EPOCH + sla_DTT (EPOCH/86400.0D0)

* Get the observatory position wrt the geocentre
      OBSLONG = ALONG(ITELNO)
      OBSLAT = ALAT(ITELNO)
      OBSHT = HGT(ITELNO)
      ST = sla_GMST (EPOCH/86400.0D0) + OBSLONG
      CALL sla_PVOBS (OBSLAT, OBSHT, ST, DPVE)

* Convert to DCS
      CALL sla_DCS2C (R2000, D2000, DPS)

* Convert UTC to TDB.
      G = (357.53 + 0.9856003 * (EPOCH / 86400.0 - 51544.5)) / RTOD
      TDB = TDT + 0.001658 * SIN (G) + 0.000014 * SIN (2 * G)
*
* Get the Earth's velocity and position using TDB
c      CALL sla_EVP (EPOCH/86400.0D0, 2000.0D0, DVB, DPB, DVH, DPH)
      CALL sla_EVP (TDB/86400.0D0, 2000.0D0, DVB, DPB, DVH, DPH)
*
* Add the EC-observatory vectors to the barycentre-EC vectors.
      DO I=1,3
        DVB(I) = DVB(I) + DPVE(I+3)
        DPB(I) = DPB(I) + DPVE(I)
      ENDDO
*
* Calculate the components of velocity and position towards pulsar.
      EVEL    = AUS*sla_DVDV(DVB,DPS)         ! sec/sec
      EDELAY  = AUS*sla_DVDV(DPB,DPS)         ! sec
      BTDB    = TDB + EDELAY
*
* Calculate the binary delay,velocity if appropriate.
      IF (PBIN.NE.0D0) THEN
        CALL EPBINARY (BTDB, EPBIN, PBIN, ASINI, WBIN, ECC,
     &                   BDELAY, BVEL, BACCN, XMA)
      ELSE
        BVEL   = 0D0
        BDELAY = 0D0
        BACCN  = 0D0
        XMA = 0D0
      ENDIF
*
* Calculate the Doppler correction
      DOPP = SQRT((1-EVEL)/(1+EVEL)*(1+BVEL)/(1-BVEL))
*
* If the period is negative, PBEPOCH, PB and PBDOT are already topocentric.
* Use these to calculate the next arrival time.
      IF (PB.LT.0D0) then
        POBS = ABS(PB)
        TOBS = ARRTIME (EPOCH,PBEPOCH,ABS(PB),PBDOT)
*
      ELSE
* Calculate the time at the pulsar and the time of the next pulse.
        PTDB  = BTDB - BDELAY
        PTIME = ARRTIME (PTDB,PBEPOCH,PB,PBDOT)
*
* PTIME-PTDB gives the first order delay of the next pulse at the
*  observatory.
        TOBS = EPOCH + (PTIME-PTDB)*DOPP
        POBS = (PB + PBDOT*(EPOCH - PBEPOCH))*DOPP
*
      ENDIF

* End of subroutine PSREPH.
      END

C
C **********************************************************************
      DOUBLE PRECISION FUNCTION ARRTIME ( EPOCH, PEPOCH, P, PDOT )
C **********************************************************************
C
C Given a reference arrival time, PEPOCH, a period, P, at that epoch, 
c     and a constant period derivative, PDOT, 
c     returns the pulse arrival time following EPOCH
c EPOCH, PEPOCH and P in secs, PDOT in sec/sec
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C Convert the periods to frequencies.
C
      V     = 1.0D0/P
      VDOT  = - PDOT / P**2
      VDDOT = 2.0D0 * PDOT**2 / P**3
C
C Calculate the number of periods between the two epochs and the
C   period at EPOCH.
C
      DT = (EPOCH - PEPOCH )
      CYCLES = V*DT + VDOT*DT**2/2.0D0 + VDDOT*(DT**3/6.0D0)
      PINST = P + PDOT*DT
C
C The offset from EPOCH to the next arrival time is the fraction 
C   of a cycle.
C
      OFFSET = MOD(CYCLES,1.0D0)
      IF ( OFFSET.LT.0.0 ) OFFSET = OFFSET + 1.0D0
C
C Calculate the arrival time.
C
      ARRTIME = EPOCH + (1.0D0-OFFSET)*PINST
C
      RETURN
C
C End of double precision function ARRTIME.
C
      END
*DECK KEPLER
C
C
C *************************************************************
      SUBROUTINE KEPLER(XMA,ECC,EA,NITS)
C *************************************************************
C
C This performs an iterative solution to keplers equation for the
C  orbit of a binary system.
C    XMA  (deg)  The mean anomaly
c    ECC         The eccentricity
C    EA   (deg)  The eccentric anomaly
C    NITS        The number of iterations used in the calculation.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EACC=1D-10,RAD=57.29577951)
C
      EA=XMA + RAD*ECC*SIND(XMA)*(1+ECC*COSD(XMA))
C
      DO 10 I=1,10
      E=(XMA + RAD*ECC*SIND(EA) - EA*ECC*COSD(EA))/(1-ECC*COSD(EA))
      IF(ABS(EA-E).LT.EACC) GO TO 20
   10 EA=E
C
   20 EA=E
      EA=(XMA + RAD*ECC*SIND(EA) - EA*ECC*COSD(EA))/(1-ECC*COSD(EA))
      NITS=I
      END

C *************************************************************
      SUBROUTINE EPBINARY (EPOCH, EPBIN, PBIN, ASINI, WBIN, ECC,
     &                       DELAY, VEL, ACCN, XMA)
C *************************************************************
*+
* Subroutine to compute the delay, velocity and acceleration of a binary
*    pulsar relative to the solar system barycentre. After Smart p358.
* All arguments are double precision.
*
* Given:
*    EPOCH       (sec)  The barycentric MJD for which data required.
*    EPBIN       (sec)  The barycentric MJD of binary periastron.
*    PBIN        (sec)  The binary barycentric period.
*    ASINI       (sec)  The binary major semi-axis.
*    WBIN        (deg)  The binary longitude of periastron.
*    ECC           -    The binary orbital eccentricity.
*
* Returned:
*    DELAY       (secs) The delay of the pulses wrt the binary b/c.
*    VEL         (s/s)  The velocity of the pulsar wrt the binary b/c.
*    ACCN        (s/s/s)The acceleration of the pulsar wrt the binary b/c.
*    XMA         (deg) Mean anomaly.
*-
*
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

*
* Define astronomical and mathematical constants.
*
      PARAMETER ( TWOPI   = 6.2831853072D0
     &           ,RTOD  = 360.0/TWOPI )
*
      E2 = (1D0 - ECC**2)
*
* Calculate the mean anomaly.
      XMA = 360.0D0/PBIN * (EPOCH - EPBIN)
*
* Solve Kepler's equation to get the eccentric anomaly.
      CALL KEPLER (XMA,ECC,E,NITS)
*
* Calculate the true anomaly.
      V = 2D0*ATAND(SQRT((1D0+ECC)/(1D0-ECC)) * TAND(E/2D0))
*
* Now calculate the required quantities
      DELAY = ASINI*(E2)*SIND(V+WBIN) / (1D0 + ECC*COSD(V))
      VEL   = TWOPI/PBIN*ASINI/SQRT(E2) * (COSD(V+WBIN)+ECC*COSD(WBIN))
      ACCN  = (TWOPI/PBIN/E2)**2*ASINI*SIND(V+WBIN)*(1D0+ECC*COSD(V))**2
*
* End of subroutine BINARY.
      END

