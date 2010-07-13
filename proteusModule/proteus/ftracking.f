C
C LAST UPDATE: 02/28/2010 (MWF)
C REMEMBER TO COMMENT OUT cc_dynamic append in gnu.py from numpy.distutils.fcompiler
C ======================================================================
C TEST F2PY FIRST
      SUBROUTINE FOO(ND,N,M,A,B)
      IMPLICIT NONE
      INTEGER ND,N,M
      DOUBLE PRECISION A(0:ND-1)
      INTEGER B(0:N-1,0:M-1)
Cf2py intent (in) ND,N,M
Cf2py double precision, intent (in) :: A(0:ND-1)
Cf2py integer intent (in) :: B(0:N-1,0:M-1)
      INTEGER I,J
      DO I=0,ND-1
         WRITE(6,*) 'A(',I,')= ',A(I)
      ENDDO
      DO I=0,N-1
         DO J=0,M-1
            WRITE(6,*) 'B(',I,',',J,')= ',B(I,J),' A(B(',I,',',J,'))= ',
     &           A(B(I,J))
         ENDDO
      ENDDO

      RETURN
      END
C TEST F2PY FIRST
      SUBROUTINE FOO2(FUN,NEQ,NNP,XG)
      IMPLICIT NONE
      EXTERNAL FUN
      INTEGER NNP,NEQ
      DOUBLE PRECISION XG(NEQ,NNP)
Cf2py intent (in) NEQ,NNP
Cf2py double precision, intent (in) :: XG(NEQ,NNP)
      INTEGER I,J
      DO I=1,NEQ
         DO J=1,NNP
            WRITE(6,*) 'XG(',I,',',J,')= ',XG(I,J)
         ENDDO
      ENDDO
      CALL FUN(NEQ)
      RETURN
      END
C
C =======================================================================
C
C MWF RELEVANT DEFINITIONS FOR PEARCE'S ORIGINAL PT123
C FORTRAN ROUTINES FOR PARTICLE TRACKING
C TO BEGIN AT LEAST, THIS WILL JUST BE REPACKAGED VERSIONS OF PEARCE'S 
C pt_adaptive_rk ROUTINES (REPACKAGED TO INTEGRATE WITH PROTEUS)
C
C < GLOBAL INFO >
C XG(I,NP) = THE I-TH COORDINATE OF THE NP-TH GLOBAL NODE
C IE(I,M) = ID OF THE GLOBAL NODE CORRESPONDING TO THE I-TH NODE OF
C           THE M-TH GLOBAL ELEMENT
C VT1E(I,J,M) = THE I-TH VELOCITY COMPONENT OF THE J-TH NODE OF THE M-TH GLOBAL ELEMENT
C               AT T1
C VT2E(I,J,M) = THE I-TH VELOCITY COMPONENT OF THE J-TH NODE OF THE M-TH GLOBAL ELEMENT
C               AT T2
C NLRL(NP) = CUMULATIVE NO. OF ELEMENTS THAT ARE CONNECTED TO GLOBAL
C            NODES 1 THROUGH NP-1
C LRL(I) = ID OF THE GLOBAL ELEMENT CORRESPONDING TO THE I-TH
C          ENTRY IN THE LRL POINT ARRAY TO PRESENT NODE-ELEMENT
C          CONNECTIVITY
C NOTE: THE NO. OF GLOBAL ELEMENTS CONNECTED TO GLOBAL NODE NP IS 
C       EQUAL TO NLRL(NP+1)-NLRL(NP).  THE ID'S OF THESE ELEMENTS ARE
C       STORED IN LRL(N1..N2), WHERE N1=NLRL(NP)+1, N2=NLRL(NP+1)
C
C < PARTICLE INFO >
C IDPT(N) = ID OF THE GLOBAL NODE CORRESPONDING TO THE N-TH TRACKED 
C           PARTICLE
C XPT(K,I,N) = THE I-TH COORDINATE OF THE K-TH LOCATION IN THE PT
C              HISTORY OF THE N-TH TRACKED PARTICLE
C
C < LOCAL INFO >
C XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF THE WORKING
C           ELEMENT
C VT1W(I,J) = THE I-TH VELOCITY COMPONENT OF THE J-TH NODE OF THE
C             WORKING ELEMENT AT T1
C VT2W(I,J) = THE I-TH VELOCITY COMPONENT OF THE J-TH NODE OF THE
C             WORKING ELEMENT AT T2
C DN_S(I) = VALUE OF THE INTERPOLATION FUNCTION ASSOCIATED WITH THE
C           I-TH NODE OF THE WORKING ELEMENT AT THE START LOCATION
C           OF A PT PROCESS
C DN(I) = VALUE OF THE INTERPOLATION FUNCTION ASSOCIATED WITH THE
C         I-TH NODE OF THE WORKING ELEMENT AT THE END LOCATION
C         OF A PT PROCESS
C
C < RK INFO >
C XS(I) = THE I-TH COORDINATE AT THE START LOCATION
C XOUT4(I) = THE I-TH COORDINATE AT THE ESTIMATED END LOCATION USING
C            THE EMBEDDED 4-TH ORDER RK
C XOUT5(I) = THE I-TH COORDINATE AT THE ESTIMATED END LOCATION USING
C            THE 5-TH ORDER RK
C XERR(I) = THE ERROR ASSOCIATED WITH THE I-TH COORDINATE WHEN XOUT4 AND
C           XOUT5 ARE COMPARED
C AK1(I), AK2(I), AK3(I), AK4(I), AK5(I), AK6(I) = VELOCITIES ESTIMATED 
C     AT VARIOUS LOCATIONS AND TIMES TO ESTIMATE XOUT4(I) AND XOUT5(I)
C XTEMP(I) = THE I-TH COORDINATE OF A WORKING LOCATION USED TO ESTIMATE 
C     THE AK FUNCTIONAL VALUES 
C ======================================================================
C
C MODEL PARAMETERS FOR SPECIFYING ARRAY SIZES
C USING ADAPTIVE RUNGE-KUTTA FOR PARTICLE TRACKING
C
C MAXNPK   = MAX. NO. OF GLOBAL NODES
C MAXELK   = MAX. NO. OF GLOBAL ELEMENTS
C MXKBDK   = MAX. NO. OF GLOBAL ELEMENTS THAT ARE CONNECTED AT A GLOBAL NODE
C MAXEQK   = MAX. NO. OF EQUATIONS TO SOLVE (1 FOR 1-D; 2 FOR 2-D; 3 FOR 3-D)
C MAXNDK   = MAX. NO. OF NODES ASSOCIATED WITH AN ELEMENT
C MAXTSK   = MAX. NO. OF TIME STEPS TAKEN INTO ACCOUNT 
C MAXPTK   = MAX. NO. OF PARTICLES POPULATED FOR TRACKING
C MAXPATHK = MAX. NO. OF TRACKING PATHS ALLOWS FOR EACH TRAKED PARTICLE
C

C ======================================================================
C MWF 
C PROTEUS TRANSLATION KEY FOR PT123
C MAXNP -- nNodes_global
C MAXEL -- nElements_global
C MXKBD -- max_nElements_node
C MAXEQ -- 3
C MAXND -- (NNDE) nNodes_element
C MAXPT -- nPointsToTrack
C MAXPATH -- not one at first, keep as max points
C            do not output path history?
C NNP  -- nNodes_global
C NEL  -- nElements_global
C NEQ  -- nSpace_global
C NPT  -- nPointsToTrack
C IBF  -- direction
C LU_O -- output file handle id
C T_START -- replaced by T_I (x_depart_times) (dimensioned for each point)
C T_END   -- replaced by TPT (x_arrive_times) (dimensioned for each point)
C DT_PT   -- replaced by T_I, TPT information
C XG      -- nodeArray,
C IE      -- elementNodesArray
C IB      -- nodeOnBoundaryArray
C NLRL    -- nodeElementOffsets
C LRL     -- nodeElementsArray
C IDPT    -- flag, a little confusing because combines node id and flag values
C XPT     -- x_out
C TPT     -- x_arrive_times
C MPT     -- element_track

      SUBROUTINE PT123
     I     (IDVE,IVERBOSE,MAXEQ,
     I      NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O,
     I      ATOL,RTOL,SF, DN_SAFE,
     I      XG,IE,NLRL,LRL,IB,
     I      VTL2G,VT1E,VT2E,T1,T2,
     I      X_I,T_I,
     M      TPT,MPT,IDPT,
     O      XPT)
     
C 
C 02/23/2010 (HPC) 
C 03/24/2010 (MWF)
C ======================================================================
C < PURPOSE > 
C   IMPLEMENT PARTICLE TRACKING USING ADAPTIVE RK ON AN ELEMENT-BY-
C   ELEMENT BASIS IN THE DESIGNATED UNSTRUTURED MESH
C MWF MODIFICATION OF PT123 FOR STEADY STATE P^1 VELOCITY
C     REPRESENTATION
C     ***NOTE***
C     EXTERNAL CODE MUST CHANGE CONNECTIVITY, LOCAL2GLOBAL, AND FLAG ARRAYS
C       TO MAKE SURE THEY ARE CONSISTENT WITH BASE ZERO 
C     THIS ROUTINE INTERNALLY SCALES THE VELOCITY BY +/- 1 BASED ON DIR
C
C
C
C
C
C TODO 
C      
C ======================================================================
C
      IMPLICIT NONE
C --HOW MUCH TO PRINT OUT 
      INTEGER IVERBOSE
Cf2py integer optional, intent (in) :: IVERBOSE = 0
C --MESH REPRESENTATION--
C THESE COULD BE PASSED IN TO MAKE DIMENSIONING SIMPLER
      INTEGER MAXEQ,MAXND
      PARAMETER(MAXND=8)
Cf2py integer optional, intent (in) :: MAXEQ = 3
C nNodes_global,nElements_global,nNodes_element,nSpace,nPointsToTrack,
C direction, output stream id
      INTEGER NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O
Cf2py integer required, intent (in) :: NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O
C NODE COORDS (nodeArray)
      DOUBLE PRECISION XG(MAXEQ*NNP)
Cf2py  double precision, intent (c), intent (in) :: XG(MAXEQ*NNP)
C ELEMENT NODE LOOKUP (elementNodesArray)
      INTEGER IE(NNDE*NEL)
Cf2py integer, intent (in) :: IE(NNDE*NEL)
C NODE - ELEMENTS IN NODE STAR LOOKUP (nodeElementOffsets,nodeElementsArray)
      INTEGER LRL(*),NLRL(*)
Cf2py integer, intent (in)  :: LRL(*),NLRL(*)
C FLAG ARRAY TO MARK NODES THAT ARE ON EXTERIOR BOUNDARY
      INTEGER IB(NNP)
Cf2py integer, intent (in)  :: IB(NNP)
C --TRACKING PARAMETERS--
      DOUBLE PRECISION ATOL,RTOL,SF,DN_SAFE
Cf2py double precision, intent (in) :: ATOL,RTOL,SF,DN_SAFE
C --TRACKING VELOCITY REPRESENTATION--
C TYPE OF LOCAL VELOCITY REPRESENTATION
C 1 -- 2 LOCAL SPACE IS C^0, P^1 
C 1 -- ASSUMED GLOBALLY CONTINUOUS (NODAL REPRESENTATION)
C 2 -- MAY BE DISCONTINUOUS (ELEMENT-BASED REPRESENTATION)
C      THIS VERSION DOES'T REALLY DISTINGUSIH SINCE 
C      VTL2G SHOULD ACCOUNT FOR THIS
C 3 -- RT0 WITH LOCAL BASIS \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}),
C 4 -- RT0 WITH LOCAL BASIS \vec N_i = \vec e_i i=0,...,d-1 and \vec N_d = \vec x
C SEE EL_VEL_PREP FOR DETAILS
      INTEGER IDVE
Cf2py integer, intent (in) :: IDVE = 2
C SEPARATE LOCAL TO GLOBAL MAP, 
C      INTEGER VTL2G(NEQ*NNDE*NEL)
CCf2py integer, intent (in) :: VTL2G(NEQ*NNDE*NEL)
      INTEGER VTL2G(*)
Cf2py integer, intent (in) :: VTL2G(*)

C DEGREES OF FREEDOM (cvelocity_dof)
      DOUBLE PRECISION VT1E(*),VT2E(*)
Cf2py double precision, intent (in) :: VT1E(*),VT2E(*)
C TIME LEVELS FOR VELOCITIES VT1E,VT2E
      DOUBLE PRECISION T1,T2
Cf2py double precision, intent (in) :: T1,T2
C

C -- INPUT TRACKING POINT INFORMATION
C POINTS TO TRACK (x_in) AND THEIR INITIAL TIMES (x_depart_times)
      DOUBLE PRECISION X_I(MAXEQ*NPT),T_I(NPT)
Cf2py double precision, intent (in) :: X_I(MAXEQ*NPT),T_I(NPT)
C
C BASE ONE, IDPT INPUT IS
C IF IDPT(K) >= 1 MEANS X(K) IS MESH NODE IDPT(K), TREAT DIFFERENTLY
C IF IDPT(K) == 0 THEN X(K) IS AN INTERIOR POINT
C IF IDPT(K) <  0 THEN DO NOT TRACK THE POINT
C ON OUTPUT SHOULD BE (BASE ONE)
C                   0  INTERIOR
C                  -1  EXITED DOMAIN SOMEWHERE 
C                  -2  DID NOT TRACK
      INTEGER IDPT(NPT)
Cf2py integer, intent (inplace) :: IDPT(NPT)
C -- IN/OUT TRACKING INFORMATION --
C TIMES TO TRACK POINTS TO (IN) FINAL TIME REACHED (OUT) 
      DOUBLE PRECISION TPT(NPT)
Cf2py double precision, intent (inplace) :: TPT(NPT)
C ELEMENT LOCATIONS FOR TRACKED POINTS
      INTEGER MPT(NPT)
Cf2py integer, intent (inplace) :: MPT(NPT)
C -- TRACKING OUTPUT
C LOCATION OF TRACKED POINTS AT END OF TRACKING CALL
      DOUBLE PRECISION XPT(MAXEQ*NPT)
Cf2py double precision, intent (inplace) :: XPT(MAXEQ*NPT)
C
C -- LOCAL VARIABLES --
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
      DOUBLE PRECISION DN_S(MAXND),DN(MAXND)
C
      DOUBLE PRECISION XOUT5(MAXEQ),XOUT4(MAXEQ),XERR(MAXEQ)
      DOUBLE PRECISION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ)
      DOUBLE PRECISION AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ)
C 
      DOUBLE PRECISION XS(MAXEQ),XTEMP(MAXEQ)
      DOUBLE PRECISION DEQ,TS,DTS,T,SDT,DT0
      DOUBLE PRECISION DIR
C     
      INTEGER I,IPT,NODE,NP,M,K
      INTEGER IDSDT,I1,I2,I3,M2,M3,N1,N2,N3,J1,J2,J3
C
C
C ===== INITIALIZATION
C
      IF (MAXEQ.NE.3) THEN
         WRITE(LU_O,*)'pt123 currently requires MAXEQ=3, but MAXEQ= ',
     &        MAXEQ
         RETURN
      ENDIF
      NODE = NNDE
      DEQ  = 1.0E0/DBLE(NEQ)
      DIR = 1.D0
      IF (IBF.EQ.-1) THEN
         DIR = -1.D0
      ENDIF
C MWF COPY X_I INTO XPT
C     HAVE TO TREAT AS FLAT ARRAY, SINCE MAY BE DIFFERENT SHAPES IN PYTHON (q, versus ebqe)
C     ASSUME MPT,TPT SET BY CALLING CODE 
CMWF DEBUG
      IF (IVERBOSE.GE.1) THEN
         WRITE(LU_O,*)'ENTERING PT123'
         WRITE(LU_O,*)'NNP= ',NNP,' NEL= ',NEL,' NNDE= ',NNDE,' NEQ= ',
     &        NEQ,' NPT= ',NPT,' IBF= ',IBF,' LU_O= ',LU_O, ' ATOL= ',
     &        ATOL,' RTOL= ',RTOL, ' SF= ',SF, 'DN_SAFE= ',DN_SAFE
      ENDIF
      DO IPT=1,NPT
         DO I=1,NEQ
            XS(I)=X_I(I + MAXEQ*(IPT-1))
            XPT(I + MAXEQ*(IPT-1))=XS(I)
        ENDDO
      ENDDO
C
C =================== START PT USING ADAPTIVE RK ====================
C
C STEP 1.  LOCATE THE VERY FIRST TIME INTERVAL FOR PT
C
C
C USE T_START AS THE REFERENCE TIME FOR PT
C
C HAVE TO SET RRTIME FOR EACH POINT NOW
C
C 1.1 FOR THE CASE OF USING A STEADY VELOCITY FIELD
C
C
C STEP 2. CONDUCT PT WITHIN THE 1-ST TIME INTERVAL
C      
      DO 600 IPT=1,NPT
        NP=IDPT(IPT)
C
C WHEN NP IS NEG, NO FURTHER TRACKING IS CONDUCTED
C
        IF(NP.LT.0)GOTO 600
C MWF NPATH NOT STORED FOR NOW
C        KPATH=NPATH(IPT)
        TS=T_I(IPT)
        DTS=TPT(IPT)-TS
        T=TS
        SDT=DTS
        DT0=SDT
        DT0=DMIN1(DT0,T2-T)
        I1=-1
        I2=-1
        I3=-1
C MWF DEBUG
        IF (IVERBOSE.GE.2) THEN
           WRITE(LU_O,*)'START OF STEP 2 IPT= ',IPT,' T1= ',T1,' T2= ',
     &          T2,' DTS= ',DTS,' TS= ',TS,' T= ',T,' NP= ',NP
        ENDIF
        DO I=1,NEQ
          XS(I)=XPT(I + MAXEQ*(IPT-1))
          IF (IVERBOSE.GE.3) THEN
C MWF DEBUG
             WRITE(LU_O,*)'XS(',I,')= ',XS(I)
          ENDIF
        ENDDO
        
        IF(NP.EQ.0)GOTO 200
C
C ===== FOR THE CASE THAT MPT(IPT)=0, I.E., THE VERY FIRST
C     TRACKING: LOOP OVER ALL CONNECTED ELEMENTS
C
        IF (IVERBOSE.GE.4) THEN
C MWF DEBUG
           WRITE(LU_O,*)' ENTERING NODE TRACKING STEP IPT= ', IPT
        ENDIF
        DO 150 I=NLRL(NP)+1,NLRL(NP+1)
          M=LRL(I)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C MWF REMOVE NPATH DIMENSIONS, CALL WITH MAXPT=NPT
          CALL ELTRAK123
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,
     I         IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I         XW,VT1W,VT2W,
     M         T,DT0,SDT,XS,
     M         AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M         XOUT5,XOUT4,XERR,DN_S,DN,
     O         IDSDT,XPT,TPT,I1,I2,I3)
C        NPATH(IPT)=KPATH
C
          IF(IDSDT.EQ.-1)THEN
            MPT(IPT)=M
C
C 0 -- INTERIOR
            IDPT(IPT)=0
            GOTO 600
          ENDIF
          IF(IDSDT.EQ.1)GOTO 250
 150    CONTINUE
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
        IF (IB(NP).EQ.1) THEN
           IDPT(IPT) = -1
        ELSE
           IDPT(IPT) = -2
        ENDIF
        IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &     (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
          WRITE(LU_O,*)'WARNING (1) IN PT123!!!' 
          WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
          WRITE(LU_O,*)'IPT= ',IPT
          WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
          WRITE(LU_O,*)'NA = NLRL(NP)+1 = ',NLRL(NP)+1
          WRITE(LU_O,*)'NB = NLRL(NP+1) = ',NLRL(NP+1)
          WRITE(LU_O,*)'LRL(NA..NB) =',(LRL(I),
     &         I=NLRL(NP)+1,NLRL(NP+1))
          WRITE(LU_O,*)
        ENDIF
        GOTO 600
C
  200   CONTINUE
C
C ===== FOR THE CASE THAT PT STARTS WITHIN AN ELEMENT,
C       I.E., MPT(IPT)=M
C
        M=MPT(IPT)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
C MWF DEBUG
        IF (IVERBOSE.GE.5) THEN
           WRITE(LU_O,*)' B4 ELTRACK ELEMENT LOOP T= ',T,' TPT= ',
     &          TPT(IPT),' IPT= ',IPT,' M= ',M
           WRITE(LU_O,*)' XS= ',(XS(I),I=1,NEQ)
        ENDIF
        CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)

C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C MWF REMOVE MAXPATH, KPATH ARGS
        CALL ELTRAK123
     I      (MAXEQ,MAXND,NPT,NEQ,NODE,
     I       IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I       XW,VT1W,VT2W,
     M       T,DT0,SDT,XS,
     M       AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M       XOUT5,XOUT4,XERR,DN_S,DN,
     O       IDSDT,XPT,TPT,I1,I2,I3)
CMWF        NPATH(IPT)=KPATH
        IF (IVERBOSE.GE.5) THEN
C MWF DEBUG
           WRITE(LU_O,*)' AFTER ELTRACK ELEMENT LOOP T= ',T,
     &          ' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT,' M= ',M
           WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
           DO K=1,NEQ
              WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',XPT(K + IPT*MAXEQ)
           ENDDO
        ENDIF
C
        IF(IDSDT.EQ.-1)THEN
          MPT(IPT)=M
C 0 -- INTERIOR
          IDPT(IPT)=0
          GOTO 600
        ENDIF
C
C ===== CONTINUE PT WITHIN A NEW ELEMENT 
C  
  250   CONTINUE
C MWF TRY TO CATCH SOMETHING WRONG WITH ELTRAK123
        IF (I1.EQ.-1.OR.I2.EQ.-1.OR.I3.EQ.-1) THEN
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
           WRITE(LU_O,*)'WARNING (5) IN PT123!!!' 
           WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
           WRITE(LU_O,*)'IPT= ',IPT
           WRITE(LU_O,*)'M= ',M
           WRITE(LU_O,*)'I1= ',I1,' I2= ',I2,' I3= ',I3
           WRITE(LU_O,*)'NP = IDPT(IPT) = ',IDPT(IPT)
C -2 -- FAILED
           IDPT(IPT) = -2
           GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON A TRIANGULAR FACE
C COMPOSED OF NODES N1, N2, N3, WHERE N1=IE(I1,M), N2=IE(I2,M), 
C N3=IE(I3,M)
C
        IF(I3.NE.0)THEN
          N1=IE(I1 + NNDE*(M-1))
          N2=IE(I2 + NNDE*(M-1))
          N3=IE(I3 + NNDE*(M-1))
          DO J1=NLRL(N1)+1,NLRL(N1+1)          
            M=LRL(J1)
            DO J2=NLRL(N2)+1,NLRL(N2+1)          
              M2=LRL(J2)
              DO J3=NLRL(N3)+1,NLRL(N3+1)          
                M3=LRL(J3)
                IF(M.EQ.M2 .AND. M2.EQ.M3)THEN
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
                   CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &                  NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &                  XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C MWF GET RID OF MAXPT,KPATH ARGS
                  CALL ELTRAK123
     I                (MAXEQ,MAXND,NPT,NEQ,NODE,
     I                 IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I                 XW,VT1W,VT2W,
     M                 T,DT0,SDT,XS,
     M                 AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M                 XOUT5,XOUT4,XERR,DN_S,DN,
     O                 IDSDT,XPT,TPT,I1,I2,I3)
C MWF                   NPATH(IPT)=KPATH
C
                  IF (IVERBOSE.GE.5) THEN
C MWF DEBUG
                     WRITE(LU_O,*)' AFTER I3.NE.0 ELTRACK ELEMENT LOOP',
     &                    'T= ',T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),
     &                    'IPT= ',IPT,' M= ',M
                     WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
                     DO K=1,NEQ
                        WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &                       XPT(K + IPT*MAXEQ)
                     ENDDO
                  ENDIF
                  IF(IDSDT.EQ.-1)THEN
                    MPT(IPT)=M
C 0 -- INTERIOR
                    IDPT(IPT)=0
                    GOTO 600
                  ENDIF
                  IF(IDSDT.EQ.1)GOTO 250
                ENDIF
              ENDDO
            ENDDO
          ENDDO
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1)THEN
            IDPT(IPT) = -1
          ELSE
            IDPT(IPT) = -2
          ENDIF
          IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &       (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
            WRITE(LU_O,*)'WARNING (2) IN PT123!!!' 
            WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
            WRITE(LU_O,*)'IPT = ',IPT
            WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
C     WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
            WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
            WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
            WRITE(LU_O,*)'N1, N2, N3 = ',N1,N2,N3
            WRITE(LU_O,*)
          ENDIF
          GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON AN EDGE COMPOSED OF 
C NODES N1 AND N2, WHERE N1=IE(I1,M), N2=IE(I2,M)
C
        IF(I2.NE.0)THEN
          N1=IE(I1 + NNDE*(M-1))
          N2=IE(I2 + NNDE*(M-1))
          DO J1=NLRL(N1)+1,NLRL(N1+1)          
            M=LRL(J1)
            DO J2=NLRL(N2)+1,NLRL(N2+1)          
              M2=LRL(J2)
              IF(M.EQ.M2)THEN
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
                 CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &                NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &                XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C REMOVE MAXPT,MAXPATH,KPATH ARGS
                CALL ELTRAK123
     I              (MAXEQ,MAXND,NPT,NEQ,NODE,
     I               IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I               XW,VT1W,VT2W,
     M               T,DT0,SDT,XS,
     M               AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M               XOUT5,XOUT4,XERR,DN_S,DN,
     O               IDSDT,XPT,TPT,I1,I2,I3)
CMWF SKIP NPATH FOR NOW
C                NPATH(IPT)=KPATH
C
                IF (IVERBOSE.GE.5) THEN
C MWF DEBUG
                   WRITE(LU_O,*)' AFTER I2.NE.0 ELTRACK ELEMENT LOOP',
     &                  'T= ',T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),
     &                  'IPT= ',IPT,' M= ',M
                   WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
                   DO K=1,NEQ
                      WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &                     XPT(K + IPT*MAXEQ)
                   ENDDO
                ENDIF
C
                IF(IDSDT.EQ.-1)THEN
                  MPT(IPT)=M
C 0 -- INTERIOR
                  IDPT(IPT)=0
                  GOTO 600
                ENDIF
                IF(IDSDT.EQ.1)GOTO 250
              ENDIF
            ENDDO
          ENDDO
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1)THEN
            IDPT(IPT) = -1
          ELSE
            IDPT(IPT) = -2
          ENDIF
          IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &       (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
            WRITE(LU_O,*)'WARNING (3) IN PT123!!!' 
            WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
            WRITE(LU_O,*)'IPT = ',IPT
            WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
CMWF          WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
            WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
            WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
            WRITE(LU_O,*)'N1, N2 = ',N1,N2
            WRITE(LU_O,*)
          ENDIF
          GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON GLOBAL NODE N1, 
C WHERE N1=IE(I1,M)
C
C === LOOP OVER ALL CONNECTED ELEMENTS
C
        N1=IE(I1 + NNDE*(M-1))
        DO J1=NLRL(N1)+1,NLRL(N1+1)          
          M=LRL(J1)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
CMWF REMOVE MAXPT, MAXPATH ARGS, KPATH
          CALL ELTRAK123
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,
     I         IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I         XW,VT1W,VT2W,
     M         T,DT0,SDT,XS,
     M         AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M         XOUT5,XOUT4,XERR,DN_S,DN,
     O         IDSDT,XPT,TPT,I1,I2,I3)
CMWF          NPATH(IPT)=KPATH
C
          IF (IVERBOSE.GE.5) THEN
C MWF DEBUG
             WRITE(LU_O,*)' AFTER I1.NE.0 ELTRACK ELEMENT LOOP T= '
     &            ,T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT
     &            ,' M= ',M
             WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
             DO K=1,NEQ
                WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',XPT(K + IPT*MAXEQ)
             ENDDO
          ENDIF
C
          IF(IDSDT.EQ.-1)THEN
            MPT(IPT)=M
C 0 -- INTERIOR
            IDPT(IPT)=0
            GOTO 600
          ENDIF
          IF(IDSDT.EQ.1)GOTO 250
        ENDDO  
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
        IF(IB(N1).EQ.1) THEN
          IDPT(IPT) = -1
        ELSE
          IDPT(IPT)=-2
        ENDIF
        IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &     (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
          WRITE(LU_O,*)'WARNING (4) IN PT123!!!' 
          WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
          WRITE(LU_O,*)'IPT = ',IPT
          WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
C        WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
          WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
          WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
          WRITE(LU_O,*)'N1 = ',N1
          WRITE(LU_O,*)
        ENDIF
C     
  600 CONTINUE
        
      RETURN
      END
C
C 
C
      SUBROUTINE ELTRAK123
     I    (MAXEQ,MAXND,MAXPT,NEQ,NODE,
     I     IPT, T1,T2, DEQ, ATOL,RTOL,SF, DN_SAFE,
     I     XW,VT1W,VT2W,
     M     T,DT0,SDT,XS,
     M     AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M     XOUT5,XOUT4,XERR,DN_S,DN,
     O     IDSDT,XPT,TPT,I1,I2,I3)
C 
C 02/08/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C IMPLEMENT PARTICLE TRACKING USING ADAPTIVE RK WITHIN AN ELEMENT
C MWF MODIFIED TO USE WITH PROTEUS
C MWF CHANGES
C     DIMENSION XPT,TPT WITHOUT NUMBER OF PATHS PER POINT  
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
CMWF NOW HAVE TO REFERENCE XPT AS FLAT ARRAY 
      DIMENSION XPT(MAXEQ*MAXPT),TPT(MAXPT)
      DIMENSION XS(MAXEQ)
      DIMENSION XOUT5(MAXEQ),XOUT4(MAXEQ),XERR(MAXEQ)
      DIMENSION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ),
     >          AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ),
     >          XTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND)
      DIMENSION ICHECK(8)
C
C =================== START PT USING ADAPTIVE RK ====================
C
C TRACKING IN ELEMENT M
C
      ICOUNT=0
      IDSDT=1
      DT=DT0
C
  100 CONTINUE
      ICOUNT=ICOUNT+1
      CALL RKCK_PT
     I    (MAXEQ,MAXND,NEQ,NODE, T,DT,T1,T2,
     I     XW,VT1W,VT2W,XS,
     M     AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     O     XOUT5,XOUT4,XERR,DN_S)
CMWF DEBUG
C      WRITE(6,*)'ELTRAK AFTER RKCK_PT T= ',T,' DT= ',DT,
C     &     ' T1= ',T1,' T2= ',T2
C      WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
C      WRITE(6,*)'DN_S= ',(DN_S(I),I=1,NODE)
C
C CHECK ERROR IN ALL THREE DIRECTIONS 
C
      RATIO=0.0D0
      DO I=1,NEQ
        XERRABS=DABS(XERR(I))
        XOUT5ABS=DABS(XOUT5(I)-XS(I))
        XOUT4ABS=DABS(XOUT4(I)-XS(I))
        XOUTMAX=DMAX1(XOUT5ABS,XOUT4ABS)
        RATIOI=XERRABS/(RTOL*XOUTMAX+ATOL)
        RATIO=DMAX1(RATIO,RATIOI)
      ENDDO
      RR=1.0E0/RATIO
C
C CASE 1: WHEN RATIO IS GREATER THAN 1
C      ==> DESIRED ACCURACY IS NOT REACHED
C      ==> TIMESTEP NEEDS TO BE REDUCED
C
      IF(RATIO.GT.1.0E0)THEN
CMWF DEBUG
C        WRITE(6,*)'ELTRAK AFTER ERROR FAILURE T= ',T,' RATIO= ',RATIO,
C     &       ' DT= ',DT,' DTT= ',DT*SF*((RR)**(0.25E0)) 
        DTT=DT*SF*((RR)**(0.25E0)) 
        DT=DTT
        GOTO 100
C
C CASE 2: WHEN RATIO IS LESS THAN OR EQUAL TO 1
C      ==> DESIRED ACCURACY IS REACHED
C      ==> TIMESTEP CAN BE INCREASED
C        
      ELSEIF(RATIO.LE.1.0E0)THEN
CMWF DEBUG
C        WRITE(6,*)'ELTRAK AFTER ERROR SUCCESS T= ',T,' RATIO= ',RATIO

C
C === EXAMINE THE COMPUTED ENDING LOCATION
C
        XSI=0.0D0
        CALL INTRP123
     I      (MAXEQ,MAXND,NEQ,NODE, XOUT5, XW,
     O       DN,DJAC)
C MWF WHAT IF NEGATIVE?
C        DL=(DJAC)**(DEQ)
        DL=DABS(DJAC)**(DEQ)
CMWF DEBUG
C        WRITE(6,*)'XOUT5= ',(XOUT5(I),I=1,NEQ)
C        WRITE(6,*)'DJAC= ',DJAC, 'DN= ',(DN(I),I=1,NODE)

C
C === ADJUST DN AND XOUT5 WHEN THE END LOCATION IS VERY CLOSE TO 
C     ELEMENT BOUNDARY
C
        JCOUNT=0
        JOUTSIDE=0
        DD=0.0E0
        DO J=1,NODE
          ICHECK(J)=0
          IF(DN(J).LT.0.0E0)THEN
            IF(DABS(DN(J)).LE.DN_SAFE)THEN
              DD=DD+DABS(DN(J))
              DN(J)=0.0E0
              ICHECK(J)=1
              JCOUNT=JCOUNT+1
            ELSE
              JOUTSIDE=JOUTSIDE+1
            ENDIF 
          ELSEIF(DN(J).GT.1.0E0)THEN
            IF(DABS(DN(J)-1.0E0).LE.DN_SAFE)THEN
              DD=DD+DABS(DN(J)-1.0E0)
              DN(J)=1.0E0
              ICHECK(J)=1
              JCOUNT=JCOUNT+1
            ELSE
              JOUTSIDE=JOUTSIDE+1
            ENDIF
          ENDIF
        ENDDO
C
        IF(JCOUNT.NE.0 .AND. JOUTSIDE.EQ.0)THEN
          DO J=1,NODE
            IF(ICHECK(J).EQ.0)THEN
              DN(J)=DN(J)-DD
              IF(DN(J).LT.0.0E0)THEN
                DD=-DN(J)
                DN(J)=0.0E0
              ELSEIF(DN(J).GT.1.0E0)THEN
                DD=DN(J)-1.0E0
                DN(J)=1.0E0
              ENDIF
            ENDIF
          ENDDO
          DO I=1,NEQ
            XOUT5(I)=0.0E0
            DO J=1,NODE
              XOUT5(I)=XOUT5(I)+DN(J)*XW(I,J)
            ENDDO
          ENDDO
        ENDIF                 
C
C === EXAMINE THE VALUES OF DN
C
CMWF DEBUG
C        WRITE(6,*)' ELTRAK EXAMINE DN, DN_S VALUES'
C        WRITE(6,*)'DN= ',(DN(II),II=1,NODE)
C        WRITE(6,*)'DN_S= ',(DN_S(II),II=1,NODE)
        DO 150 I=1,NODE

C
C A. WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C    ===> COMPUTE XSI
C
          IF(DN(I).LT.0.0E0)THEN
C MWF DEBUG
C            WRITE(6,*)' DN(',I,')= ',DN(I),' OUTSIDE ELEMENT'

            D1=DN_S(I)
            D2=DN(I)
            D12=DN_S(I)-DN(I)
C
C A1. WHEN THERE IS NO TRACKIING THROUGH THIS ELEMENT
C
            IF(DABS(D1*DL).LT.ATOL)THEN
C
C IF THE PT IS LEAVING THE ELEMENT
C ==> SET IDSDT TO 0, AND IDENTIFY I1, I2, I3 TO
C     CONDUCT PT IN AN ADJACENT ELEMENT
C
              IF(DABS(D2*DL).GT.ATOL)THEN
                CALL DN_CHECK
     I              (MAXND,NODE,NEQ, DL,ATOL,
     I               DN_S,
     O               I1,I2,I3)
                IDSDT=0
C MWF DEBUG
C                WRITE(6,*)' ELTRAK DN_CHECK LEAVING IDSDT= ',IDSDT, 
C     &               ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
C     &               ' D1= ',D1, ' D2= ',D2,' D12= ',D12
                RETURN
              ENDIF
C
C A2. WHEN THERE IS A TRACKING THROUGH THIS ELEMENT
C
            ELSE
C
C IF THE ENDING LOCATION CAN BE CONSIDERED ON THE ELEMENT BOUNDARY
C ==> NO NEED TO REDUCE TIMESTEP
C
              IF(DABS(D2*DL).LE.ATOL)GOTO 150
C
C IF THE ENDING LOCATION IS TRULY OUTSIDE OF THE ELEMENT
C ==> COMPUTE AN INTERPOLATION FACTOR
C
              XSI=DMAX1(XSI,D12/D1)
            ENDIF
          ENDIF
  150   CONTINUE
C
C === IF XSI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C
        IF(XSI.GT.1.0E0)THEN
          DTT=DT/XSI
          DT=DTT
          GOTO 100
        ENDIF
C
C B. WHEN THE ENDING LOCATION IS EITHER WITHIN THE ELEMENT OR 
C    ON THE ELEMENT BOUDNARY
C ==> UPDATE INFORMATION FOR THE SUCCESSIVE PT
C
        T=T+DT
CMWF        KPATH=KPATH+1
CMWF        TPT(KPATH,IPT)=T
        TPT(IPT)=T
        DO I=1,NEQ
          XS(I)=XOUT5(I)
CMWF          XPT(KPATH,I,IPT)=XOUT5(I)
CMWF NOW FLATTEN INDEXING
          XPT(I + MAXEQ*(IPT-1))=XOUT5(I)
        ENDDO   
        DTT=DT*SF*((RR)**(0.2E0))
        SDT=SDT-DT
        DT0=DT0-DT
C
C IF THE TRACKING TIME IS COMPLETELY CONSUMED
C ==> SET IDSDT TO -1
C
        IF(DABS(SDT).LE.1.0E-10)THEN
          IDSDT=-1
          SDT=0.0E0
        ENDIF
        IF(DABS(DT0).LE.1.0E-10)THEN
          IDSDT=-1
          DT0=0.0E0
        ENDIF
C
C IF THE ENDING LOCATION IS ON THE ELEMENT BOUNDARY
C ==> EXIT PT IN THIS ELEMENT
C
        DO I=1,NODE
          IF(DABS(DN(I)*DL).LT.ATOL)THEN
            CALL DN_CHECK
     I          (MAXND,NODE,NEQ, DL,ATOL,
     I           DN,
     O           I1,I2,I3)
C MWF DEBUG
C            WRITE(6,*)' ELTRAK DN_CHECK ON BNDY STAYING IDSDT= ',IDSDT, 
C     &           ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
C     &           ' DN_S= ',(DN_S(II),II=1,NODE)
            RETURN
          ENDIF
        ENDDO
C
C UPDATE THE TRACKING TIMESTEP (I.E., DT) FOR THE SUCCESSIVE TRACKING
C WITHIN THIS SAME ELEMENT
C
CMWF ORIG
C        IF(IDSDT.EQ.-1)RETURN
C MWF DEBUG
        IF(IDSDT.EQ.-1) THEN
C MWF DEBUG
C            WRITE(6,*)' ELTRAK DN_CHECK ON BNDY STAYING IDSDT= ',IDSDT, 
C     &           ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
C     &           ' DN_S= ',(DN_S(II),II=1,NODE)
           RETURN
        ENDIF
C MWF END DEBUG
        DT=DT0
        DT=DMIN1(DT,DTT)
        ICOUNT=0
        GOTO 100
C
      ENDIF
C 
C  999 CONTINUE
      RETURN
      END

C
C 
C
      SUBROUTINE DN_CHECK
     I    (MAXND,NODE,NEQ, DL,ATOL,
     I     DN,
     O     I1,I2,I3)
C 
C 02/03/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C PREPARE I1, I2, I3 FOR THE SUCCESSIVE PT WHEN THE TRACKED NODE IS
C ON THE ELEMENT BOUNDARY
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DN(MAXND)
C
C === CHECK IF THE PT ENDS ON ANY SIDE OF THE ELEMENT
C
      I1=0
      I2=0
      I3=0
      DO I=1,NODE
        DD=DN(I)*DL
        IF(DABS(DD).LE.ATOL)THEN
          IF(I1.EQ.0)THEN
            I1=I
          ELSE
            IF(I2.EQ.0)THEN
              I2=I
            ELSE
              I3=I
            ENDIF
          ENDIF
        ENDIF                
      ENDDO
C
C === IF THREE INTERPOLATION FUNCTIONS ARE COMPUTED TO BE
C     ZERO 
C
C THIS WILL OCCUR ONLY IN 3-D: THE TRACKED NODE COINCIDES WITH
C                              AN ELEMENT NODE (I.E., I1)
C
      IF(I3.NE.0)THEN
        DO I=1,NODE
          IF(I.NE.I1 .AND. I.NE.I2 .AND. I.NE.I3)THEN
            I1=I
            I2=0
            I3=0
            GOTO 999
          ENDIF
        ENDDO         
      ENDIF
C
C === IF TWO INTERPOLATION FUNCTIONS ARE COMPUTED TO BE
C     ZERO 
C
      IF(I2.NE.0)THEN
C
C FOR THE CASE OF 2-D: THE TRACKED NODE COINCIDES WITH
C                      AN ELEMENT NODE (I.E., I1)
C
        IF(NEQ.EQ.2)THEN
          DO I=1,NODE
            IF(I.NE.I1 .AND. I.NE.I2)THEN
              I1=I
              I2=0
              I3=0
              GOTO 999
            ENDIF
          ENDDO
C
C FOR THE CASE OF 3-D: THE TRACKED NODE FALL ON AN
C                      ELEMENT EDGE (I.E., I1-I2)
C
        ELSEIF(NEQ.EQ.3)THEN
          DO I=1,NODE
            IF(I.NE.I1 .AND. I.NE.I2)THEN           
              IF(I3.EQ.0)THEN
                I3=I
              ELSE
                I1=I3
                I2=I
                I3=0
                GOTO 999
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C === IF ONE INTERPOLATION FUNCTIONS ARE COMPUTED TO BE
C     ZERO 
C
      IF(I1.NE.0)THEN
C
C FOR THE CASE OF 1-D: THE TRACKED NODE COINCIDES WITH
C                      AN ELEMENT NODE (I.E., I1)
C
        IF(NEQ.EQ.1)THEN
          DO I=1,NODE
            IF(I.NE.I1)THEN
              I1=I
              I2=0
              I3=0
              GOTO 999
            ENDIF
          ENDDO
C
C FOR THE CASE OF 2-D: THE TRACKED NODE FALL ON AN
C                      ELEMENT EDGE (I.E., I1-I2)
C
        ELSEIF(NEQ.EQ.2)THEN
          DO I=1,NODE
            IF(I.NE.I1)THEN           
              IF(I3.EQ.0)THEN
                I3=I
              ELSE
                I1=I3
                I2=I
                I3=0
                GOTO 999
              ENDIF
            ENDIF
          ENDDO
C
C FOR THE CASE OF 3-D: THE TRACKED NODE FALL ON AN
C                      ELEMENT FACE (I.E., I1-I2-I3)
C
        ELSEIF(NEQ.EQ.3)THEN
          DO I=1,NODE
            IF(I.NE.I1)THEN           
              IF(I3.EQ.0)THEN
                I3=I
              ELSE
                IF(I2.EQ.0)THEN
                  I2=I
                ELSE
                  I1=I
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF           
      ENDIF  
C
  999 CONTINUE
      RETURN
      END 

      SUBROUTINE RKCK_PT
     I    (MAXEQ,MAXND,NEQ,NODE, T,DT, T1,T2,
     I     XW,VT1W,VT2W, XS,
     M     AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     O     XOUT5,XOUT4,XERR,DN_S)
C 
C 02/08/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C   ESTIMATE ERRORS BETWEEN 4TH- AND 5TH-ORDER RUNGE-KUTTA (RK) USING
C   CASH-KARP RK METHOD CONCERNING PARTICLE TRACKING (PK)
C < INPUT > 
C   MAXEQ  = MAX. NO. OF 1ST-ORDER ODE
C   NEQ    = NO. OF 1ST-ORDER ODE
C   T      = STARTING TIME FOR PK
C   DT     = AVAILABLE TRACKING TIME
C   XS     = START LOCATION 
C   T1, T2 = TIMES WHERE VELOCITIES (TIME DERIVATIVES) ARE GIVEN
C   XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF A SPECIFIC 
C             ELEMENT
C   VT1W(I,J) = THE I-TH VELOCITY COMPONENT ASSOCAITED WITH THE J-TH 
C               NODE AT TIME T1
C   VT2W(I,J) = THE I-TH VELOCITY COMPONENT ASSOCAITED WITH THE J-TH 
C               NODE AT TIME T2
C < OUTPUT >
C   XOUT5  = ESTIMATE FROM 5TH-ORDER RK
C   XOUT4  = ESTIMATE FROM 4TH-ORDER RK
C   XERR   = ERROR ESTIMATE BETWEEN 4TH- AND 5TH-ORDER RK
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XS(MAXEQ)
      DIMENSION XOUT5(MAXEQ),XOUT4(MAXEQ),XERR(MAXEQ)
      DIMENSION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ),
     >          AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ),
     >          XTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND)
C
C ===== DEFINE THE PARAMETERS USED FOR CASH-KARP RK
C
      A2  = 0.2E0
      A3  = 0.3E0
      A4  = 0.6E0
      A5  = 1.0E0
      A6  = 0.875E0
C
      B21 = 0.2E0
      B31 = 0.075E0
      B41 = 0.3E0
      B51 = -11.E0/54.E0
      B61 = 1631.E0/55296.E0
      B32 = 0.225E0
      B42 = -0.9E0
      B52 = 2.5E0
      B62 = 175.E0/512.E0
      B43 = 1.2E0
      B53 = -70.E0/27.E0
      B63 = 575.E0/13824.E0
      B54 = 35.E0/27.E0
      B64 = 44275.E0/110592.E0
      B65 = 253.E0/4096.E0
C
      C1  = 37.E0/378.E0
      C2  = 0.0E0
      C3  = 250.E0/621.E0
      C4  = 125.E0/594.E0
      C5  = 0.0E0
      C6  = 512.E0/1771.E0
C
      D1  = 2825.E0/27648.E0
      D2  = 0.0E0
      D3  = 18575.E0/48384.E0
      D4  = 13525.E0/55296.E0
      D5  = 277.E0/14336.E0
      D6  = 0.25E0
C
      E1  = C1-D1
      E2  = C2-D2
      E3  = C3-D3
      E4  = C4-D4
      E5  = C5-D5
      E6  = C6-D6
C
C INITIALIZE AK2, AK3, AK4, AK5, AK6
C
      DO I=1,NEQ
        AK2(I)=0.0E0
        AK3(I)=0.0E0
        AK4(I)=0.0E0
        AK5(I)=0.0E0
        AK6(I)=0.0E0
      ENDDO
C
C COMPUTE TIME DERIVATIVE FUNCTIONAL VALUES NEEDED FOR CASH-KARP RK
C THESE FUNCTIONS ARE AK2, AK3, AK4, AK5, AND AK6
C
C STEP 0:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)
      ENDDO
      TT=T
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK1,DN_S)
C
C STEP 1:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*B21*AK1(I)
      ENDDO
      TT=T+A2*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK2,DN)
C
C STEP 2:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B31*AK1(I)+B32*AK2(I))
      ENDDO
      TT=T+A3*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK3,DN)
C
C STEP 3:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B41*AK1(I)+B42*AK2(I)+B43*AK3(I))
      ENDDO
      TT=T+A4*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK4,DN)
C
C STEP 4:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B51*AK1(I)+B52*AK2(I)+B53*AK3(I)+
     >                     B54*AK4(I))
      ENDDO
      TT=T+A5*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK5,DN)
C
C STEP 5:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B61*AK1(I)+B62*AK2(I)+B63*AK3(I)+
     >                     B64*AK4(I)+B65*AK5(I))
      ENDDO
      TT=T+A6*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK6,DN)
C
C ESTIMATE ERROR USING 5TH- AND 4TH-ORDER RK
C
      DO I=1,NEQ
        XOUT5(I)=XS(I)+DT*(C1*AK1(I)+C2*AK2(I)+C3*AK3(I)+
     >                     C4*AK4(I)+C5*AK5(I)+C6*AK6(I))
        XOUT4(I)=XS(I)+DT*(D1*AK1(I)+D2*AK2(I)+D3*AK3(I)+
     >                     D4*AK4(I)+D5*AK5(I)+D6*AK6(I))
        XERR(I)=DT*(E1*AK1(I)+E2*AK2(I)+E3*AK3(I)+
     >              E4*AK4(I)+E5*AK5(I)+E6*AK6(I))
      ENDDO
C
C ===== RETURN TO THE CALLING ROUTINE
C
      RETURN
      END
C
C
C
      SUBROUTINE VEL123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W,
     O     AK,DN)
C 
C 02/08/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C   COMPUTE TIME DERIVATIVE FUNCTIONAL VALUES NEEDED FOR CASH-KARP RK
C < INPUT > 
C   MAXEQ  = MAX. NO. OF 1ST-ORDER ODE
C   NEQ     = NO. OF 1ST-ORDER ODE
C   TT      = TIME USED FOR COMPUTATION
C   XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF A SPECIFIC 
C             ELEMENT
C   VT1W(I,J) = THE I-TH VELOCITY COMPONENT ASSOCAITED WITH THE J-TH 
C               NODE AT TIME T1
C   VT2W(I,J) = THE I-TH VELOCITY COMPONENT ASSOCAITED WITH THE J-TH 
C               NODE AT TIME T2
C < OUTPUT >
C   AK      = FUNCTIONAL VALUES OUT OF COMPUTATION
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),AK(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN(MAXND)
      DIMENSION VT1(3),VT2(3)
C
C
C === COMPUTE THE INTERPOLATION FUNCTIONAL VALUES
C
C IN SPACE:
C
      CALL INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP, XW,
     O     DN,DJAC)
C
C IN TIME:
C
CMWF ALLOW T2=T1?
      IF(DABS(T2-T1).LE.1.0D-9) THEN
         ETA=0.D0
      ELSE
         ETA=(TT-T1)/(T2-T1)
      ENDIF
C
C
C === DO INTERPOLATION
C
      DO I=1,NEQ
        VT1(I)=0.0E0
        VT2(I)=0.0E0
      ENDDO
      DO I=1,NEQ
        DO J=1,NODE
          VT1(I)=VT1(I)+VT1W(I,J)*DN(J)
          VT2(I)=VT2(I)+VT2W(I,J)*DN(J)
        ENDDO
      ENDDO
      DO I=1,NEQ
        AK(I)=VT1(I)+ETA*(VT2(I)-VT1(I))
      ENDDO    
C MWF DEBUG
      WRITE(6,*)' VEL123 T1= ',T1, 'T2= ',T2,' TT= ',TT,' ETA= ',ETA
      DO I=1,NEQ
         WRITE(6,*)' VT1(',I,')= ',VT1(I),' VT2(',I,')= ',VT2(I),' AK(',
     &        I,')= ',AK(I)
      ENDDO
C
C  999 CONTINUE
      RETURN
      END 
C
C
C
      SUBROUTINE INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP, XW,
     O     DN,DJAC)
C
C 02/01/2010 (HPC) 
C ======================================================================
C < PURPOSE > 
C   COMPUTE THE VALUES OF INTERPOLATION FUNCTIONS
C < INPUT > 
C   XTEMP(I) = THE I-TH COORDINATE THE THE LOCATION USED FOR COMPUTATION
C   XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF A SPECIFIC 
C             ELEMENT
C < OUTPUT >
C   DN(I) = VALUE ASSCIATED WITH THE INTERPOLATION FUNCTION ASSOCIATED
C           WITH THE I-TH NODE
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
      DIMENSION DN(MAXND)
      DIMENSION A(4),B(4),C(4),D(4)
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
        DJAC=XW(1,2)-XW(1,1)
        DN(2)=(XTEMP(1)-XW(1,1))/DJAC
        DN(1)=1.0E0-DN(2)
C
C
C ===== FOR THE CASE OF A 2-D TRIANGULAR ELEMENT
C
      ELSEIF(NEQ.EQ.2)THEN
        X12=XW(1,1)-XW(1,2)
        X23=XW(1,2)-XW(1,3)
        X31=XW(1,3)-XW(1,1)  
        Y12=XW(2,1)-XW(2,2)
        Y23=XW(2,2)-XW(2,3)
        Y31=XW(2,3)-XW(2,1)
        DJAC=XW(1,1)*Y23+XW(1,2)*Y31+XW(1,3)*Y12
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
C        DJAC=DABS(DJAC)
        DN(1)=Y23*XTEMP(1)-X23*XTEMP(2)+XW(1,2)*XW(2,3)-XW(1,3)*XW(2,2)
        DN(2)=Y31*XTEMP(1)-X31*XTEMP(2)+XW(1,3)*XW(2,1)-XW(1,1)*XW(2,3)
        DN(3)=Y12*XTEMP(1)-X12*XTEMP(2)+XW(1,1)*XW(2,2)-XW(1,2)*XW(2,1)
        DO I=1,NODE
          DN(I)=DN(I)/DJAC
        ENDDO
C
C
C ===== FOR THE CASE OF A 3-D TETRAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3)THEN
        DJAC=0.0D0
        DO KK=1,4
          IF(KK.EQ.1)THEN
            K1=2
            K2=3
            K3=4
          ELSEIF(KK.EQ.2)THEN
            K1=1
            K2=3
            K3=4
          ELSEIF(KK.EQ.3)THEN
            K1=1
            K2=2
            K3=4
          ELSE
            K1=1
            K2=2
            K3=3
          ENDIF
C
          A(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(2,K2)*XW(3,K3)+
     >          XW(2,K1)*XW(3,K2)*XW(1,K3)+XW(3,K1)*XW(1,K2)*XW(2,K3)-
     >          XW(1,K3)*XW(2,K2)*XW(3,K1)-XW(2,K3)*XW(3,K2)*XW(1,K1)-
     >          XW(3,K3)*XW(1,K2)*XW(2,K1))
C
          B(KK)=(-1.0D0)**KK*(XW(2,K1)*XW(3,K2)+
     >          XW(2,K2)*XW(3,K3)+XW(2,K3)*XW(3,K1)-
     >          XW(2,K3)*XW(3,K2)-XW(2,K2)*XW(3,K1)-
     >          XW(2,K1)*XW(3,K3))
C
          C(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(3,K2)+
     >          XW(1,K2)*XW(3,K3)+XW(1,K3)*XW(3,K1)-
     >          XW(1,K3)*XW(3,K2)-XW(1,K2)*XW(3,K1)-
     >          XW(1,K1)*XW(3,K3))
C
          D(KK)=(-1.0D0)**KK*(XW(1,K1)*XW(2,K2)+
     >          XW(1,K2)*XW(2,K3)+XW(1,K3)*XW(2,K1)-
     >          XW(1,K3)*XW(2,K2)-XW(1,K2)*XW(2,K1)-
     >          XW(1,K1)*XW(2,K3))
C
          DJAC=DJAC+A(KK)
          DN(KK)=A(KK)+B(KK)*XTEMP(1)+C(KK)*XTEMP(2)+D(KK)*XTEMP(3)
        ENDDO
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
C        DJAC=DABS(DJAC)
        DO KK=1,4
          DN(KK)=DN(KK)/DJAC
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
C

      SUBROUTINE EL_VEL_PREP
     I    (MAXND,MAXEQ,NNDE,
     I     NNP,NEL,NODE,NEQ, M,IDVE,DIR,
     I     XG,IE,VTL2G,VT1E,VT2E,
     O     XW,VT1W,VT2W)
C 
C 02/23/2010 (HPC)
C 05/10/2010 (MWF)
C ======================================================================
C < PURPOSE > 
C PREPARE ELEMENT NODAL VELOCITY FOR PT WITHIN THE ELEMENT
C MWF MODIFIED PEARCE'S ORIGINAL ROUTINE TO TAKE A LOCAL TO GLOBAL DOF
C     MAP AND ALLOW RT0 (ON SIMPLICES) AS WELL
C ======================================================================
C
C      IMPLICIT REAL*8(A-H,O-Z)
C
C      DIMENSION XG(MAXEQ,MAXNP),IE(MAXND,MAXEL)
C      DIMENSION VT1N(MAXEQ,MAXNP),VT2N(MAXEQ,MAXNP)
C      DIMENSION VT1E(MAXEQ,MAXND,MAXEL),VT2E(MAXEQ,MAXND,MAXEL)
C      DIMENSION XW(MAXEQ,MAXND)
C      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
C  
      IMPLICIT NONE
C DIMENSIONS FOR MAXIMUM 
C  MAXIMUM NUMBER OF NODES PER ELEMENT (8)
C  MAXIMUM NUMBER OF EQUATIONS FOR VELOCITY (3)   
      INTEGER MAXND,MAXEQ
C ACTUAL DIMENSIONS FOR NUMBER OF NODES, NUMBER OF ELEMENTS, 
C NUMBER OF NODES PER ELEMENT (MAX FOR THIS MESH ...) 
C NUMBER OF NODES FOR THIS ELEMENT, AND NUMBER OF EQUATIONS TO 
C INTEGRATE FOR VELOCITY (SPACE DIMENSION)
C THIS ELEMENT IF HAVE MIXED TYPES
      INTEGER NNP,NEL,NNDE,NODE,NEQ
C CURRENT ELEMENT
      INTEGER M
C FLAG FOR VELOCITY TYPE
C 1 -- 2 LOCAL SPACE IS C^0, P^1 
C 1 -- ASSUMED GLOBALLY CONTINUOUS (NODAL REPRESENTATION)
C 2 -- MAY BE DISCONTINUOUS (ELEMENT-BASED REPRESENTATION)
C      THIS VERSION DOES'T REALLY DISTINGUSIH SINCE 
C      VTL2G SHOULD ACCOUNT FOR THIS
C 3 -- RT0 WITH LOCAL BASIS \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}),
C 4 -- RT0 WITH LOCAL BASIS \vec N_i = \vec e_i i=0,...,d-1 and \vec N_d = \vec x
      INTEGER IDVE
C DIRECTION IN TIME FOR SCALING VELOCITY
      DOUBLE PRECISION DIR
 
C NODAL COORDINATES
      DOUBLE PRECISION XG(MAXEQ,NNP)
C ELEMENT -- NODES CONNECTIVITY
      INTEGER IE(NNDE,NEL)
C LOCAL TO GLOBAL DOF MAPPING
C      INTEGER VTL2G(NEQ,NNDE,NEL)
      INTEGER VTL2G(*)
C DEGREES OF FREEDOM AT TIME LEVELS 1 AND 2
      DOUBLE PRECISION VT1E(*),VT2E(*)
C WORK ARRAY FOR NODES AND VELOCITIES AT TIME LEVEL 1 AND 2
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)

C LOCAL VARIABLES
      INTEGER I,J,K,IVDOF,IEM
      DOUBLE PRECISION DVOL,DRT0FACT
C
C 
C WHEN IDVE = 1,2 USE C0P1 VELOCITY 
C
      IF(IDVE.EQ.1.OR.IDVE.EQ.2)THEN
        DO J=1,NODE
          IEM = IE(J,M)
          DO K=1,NEQ
            XW(K,J)=XG(K,IEM)
C            IVDOF = VTL2G(K,J,M)
            IVDOF = VTL2G(K + NEQ*(J-1) + NEQ*NNDE*(M-1))
            VT1W(K,J)=VT1E(IVDOF)*DIR
            VT2W(K,J)=VT2E(IVDOF)*DIR
C MWF DEBUG
            WRITE(6,*)'EL_VEL_PREP IDVE.EQ.2 M= ',M,' J= ',J,
     &           ' IEM= ',IEM,' XW(',K,',',J,')= ',XW(K,J),
     &           ' IVDOF= ',IVDOF,' VT1W(',K,',',J,')= ',VT1W(K,J),
     &           ' VT2W(',K,',',J,')= ',VT2W(K,J)
C MWF DEBUG            
          ENDDO
        ENDDO
      ELSE IF (IDVE.EQ.3) THEN
C RT0 VELOCITY ON A SIMPLEX WITH REPRESENTATION
C \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}), i=0,...,d
C  d IS THE SPACE DIMENSION, |\Omege_e| IS THE ELEMENT VOLUMNE
C AND \vec x_{n,i} IS THE NODE ACROSS FROM FACE i
C LOAD LOCAL NODAL COORDS FIRST
        DO J=1,NODE
          IEM = IE(J,M)
          DO K=1,NEQ
            XW(K,J)=XG(K,IEM)
          ENDDO
        ENDDO
        CALL ELEMENT_VOLUME(MAXEQ,MAXND,NEQ,NODE,XW,DVOL)
        DRT0FACT = 1.D0/(DBLE(NEQ)*DVOL)
        DO J=1,NODE
          DO K=1,NEQ
            VT1W(K,J)=0.D0
            VT2W(K,J)=0.D0
            DO I=1,NODE
C TODO TRY TO USE A SUBSET OF VTL2G MEMORY FOR VECTOR VALUED SPACES
C      OR JUST MAKE SIZING CONSISTENT WITH NLOCAL_DOF X NELEMENTS?
              IVDOF = VTL2G(I + NNDE*(M-1))
              VT1W(K,J) = VT1W(K,J) + VT1E(IVDOF)*
     &             DRT0FACT*(XW(K,J)-XW(K,I))*DIR
              VT2W(K,J) = VT2W(K,J) + VT2E(IVDOF)*
     &             DRT0FACT*(XW(K,J)-XW(K,I))*DIR
            ENDDO
          ENDDO
        ENDDO     
      ELSE IF (IDVE.EQ.4) THEN
C RT0 VELOCITY ON A SIMPLEX WITH REPRESENTATION
C \vec v_e = \vec a_e + b_e \vec x
C WITH BASIS NUMBERING CONVENTION
C \vec N_i = \vec e_i i=0,...,d-1 and \vec N_d = \vec x
C WHERE d IS THE SPACE DIMENSION, e IS THE LOCAL ELEMENT
C LOAD LOCAL COORDINATES FIRST
        DO J=1,NODE
          IEM = IE(J,M)
          DO K=1,NEQ
            XW(K,J)=XG(K,IEM)
          ENDDO
        ENDDO
        DO J=1,NODE
          DO K=1,NEQ
            VT1W(K,J)=0.D0
            VT2W(K,J)=0.D0
            DO I=1,NODE-1
C TODO TRY TO USE A SUBSET OF VTL2G MEMORY FOR VECTOR VALUED SPACES
C      OR JUST MAKE SIZING CONSISTENT WITH NLOCAL_DOF X NELEMENTS?
C              IVDOF = VTL2G(1,I,M)
              IVDOF = VTL2G(I + NNDE*(M-1))
              VT1W(K,J)=VT1W(K,J) + VT1E(IVDOF)*DIR
              VT2W(K,J)=VT2W(K,J) + VT2E(IVDOF)*DIR
            ENDDO
            I = NODE
C            IVDOF = VTL2G(1,I,M)
            IVDOF = VTL2G(I + NNDE*(M-1))
            VT1W(K,J)=VT1W(K,J) + VT1E(IVDOF)*XW(K,J)*DIR
            VT2W(K,J)=VT2W(K,J) + VT2E(IVDOF)*XW(K,J)*DIR
          ENDDO
        ENDDO
      ELSE
        WRITE(6,*)'EL_VEL_PREP IDVE= ',IDVE,'NOT VALID QUITTING'
        CALL EXIT(1)
      ENDIF
C

      RETURN
      END
C
C
C
      SUBROUTINE ELEMENT_VOLUME
     I    (MAXEQ,MAXND,NEQ,NODE,XW,
     O     DVOL)
C
C 05/10/2010 (MWF) 
C ======================================================================
C < PURPOSE > 
C   COMPUTE VOLUME OF ELEMENT
C < INPUT > 
C   XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF A SPECIFIC 
C             ELEMENT
C < OUTPUT >
C   DVOL(I) = VOLUME OF CURRENT ELEMENT
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION A(4),B(4),C(4),D(4)
      DOUBLE PRECISION DVOL,DJAC
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
        DJAC=XW(1,2)-XW(1,1)
        DVOL=DABS(DJAC)
C
C
C ===== FOR THE CASE OF A 2-D TRIANGULAR ELEMENT
C
      ELSEIF(NEQ.EQ.2)THEN
        X12=XW(1,1)-XW(1,2)
        X23=XW(1,2)-XW(1,3)
        X31=XW(1,3)-XW(1,1)  
        Y12=XW(2,1)-XW(2,2)
        Y23=XW(2,2)-XW(2,3)
        Y31=XW(2,3)-XW(2,1)
        DJAC=XW(1,1)*Y23+XW(1,2)*Y31+XW(1,3)*Y12
        DVOL=0.5D0*DABS(DJAC)
C
C
C ===== FOR THE CASE OF A 3-D TETRAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3)THEN
        DJAC=0.0D0
        DO KK=1,4
          IF(KK.EQ.1)THEN
            K1=2
            K2=3
            K3=4
          ELSEIF(KK.EQ.2)THEN
            K1=1
            K2=3
            K3=4
          ELSEIF(KK.EQ.3)THEN
            K1=1
            K2=2
            K3=4
          ELSE
            K1=1
            K2=2
            K3=3
          ENDIF
C
          A(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(2,K2)*XW(3,K3)+
     >          XW(2,K1)*XW(3,K2)*XW(1,K3)+XW(3,K1)*XW(1,K2)*XW(2,K3)-
     >          XW(1,K3)*XW(2,K2)*XW(3,K1)-XW(2,K3)*XW(3,K2)*XW(1,K1)-
     >          XW(3,K3)*XW(1,K2)*XW(2,K1))
C
          B(KK)=(-1.0D0)**KK*(XW(2,K1)*XW(3,K2)+
     >          XW(2,K2)*XW(3,K3)+XW(2,K3)*XW(3,K1)-
     >          XW(2,K3)*XW(3,K2)-XW(2,K2)*XW(3,K1)-
     >          XW(2,K1)*XW(3,K3))
C
          C(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(3,K2)+
     >          XW(1,K2)*XW(3,K3)+XW(1,K3)*XW(3,K1)-
     >          XW(1,K3)*XW(3,K2)-XW(1,K2)*XW(3,K1)-
     >          XW(1,K1)*XW(3,K3))
C
          D(KK)=(-1.0D0)**KK*(XW(1,K1)*XW(2,K2)+
     >          XW(1,K2)*XW(2,K3)+XW(1,K3)*XW(2,K1)-
     >          XW(1,K3)*XW(2,K2)-XW(1,K2)*XW(2,K1)-
     >          XW(1,K1)*XW(2,K3))
C
          DJAC=DJAC+A(KK)
        ENDDO
        DVOL=DABS(DJAC)/6.D0
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE DERIV123
     I    (MAXEQ,MAXND,NEQ,NODE, XTEMP, XW,
     O     DDN,DJAC)
C
C 05/12/2010 (MWF) 
C ======================================================================
C < PURPOSE > 
C   COMPUTE THE DERIVATIVES OF INTERPOLATION FUNCTIONS
C    COULD PUT THIS IN INTRP123 AFTER WE ARE SURE ANALYTICAL TRACKING IS USEFUL
C < INPUT > 
C   XTEMP(I) = THE I-TH COORDINATE THE THE LOCATION USED FOR COMPUTATION
C   XW(I,J) = THE I-TH COORDINATE OF THE J-TH NODE OF A SPECIFIC 
C             ELEMENT
C < OUTPUT >
C   DDN(J,I) = DERIVATIVE WITH RESPECT TO JTH COORDINATE OF THE INTERPOLATION FUNCTION ASSOCIATED
C     WITH THE I-TH NODE
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
      DIMENSION DDN(MAXEQ,MAXND)
      DIMENSION A(4),B(4),C(4),D(4)
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
        DJAC=XW(1,2)-XW(1,1)
C        DN(2)=(XTEMP(1)-XW(1,1))/DJAC
        DDN(1,2)=1.D0/DJAC
C        DN(1)=1.0E0-DN(2)
        DDN(1,1)=-DDN(1,2)
C
C
C ===== FOR THE CASE OF A 2-D TRIANGULAR ELEMENT
C
      ELSEIF(NEQ.EQ.2)THEN
        X12=XW(1,1)-XW(1,2)
        X23=XW(1,2)-XW(1,3)
        X31=XW(1,3)-XW(1,1)  
        Y12=XW(2,1)-XW(2,2)
        Y23=XW(2,2)-XW(2,3)
        Y31=XW(2,3)-XW(2,1)
        DJAC=XW(1,1)*Y23+XW(1,2)*Y31+XW(1,3)*Y12
C        DN(1)=Y23*XTEMP(1)-X23*XTEMP(2)+XW(1,2)*XW(2,3)-XW(1,3)*XW(2,2)
        DDN(1,1)= Y23
        DDN(2,1)=-X23
C        DN(2)=Y31*XTEMP(1)-X31*XTEMP(2)+XW(1,3)*XW(2,1)-XW(1,1)*XW(2,3)
        DDN(1,2)= Y31
        DDN(2,2)=-X31
C        DN(3)=Y12*XTEMP(1)-X12*XTEMP(2)+XW(1,1)*XW(2,2)-XW(1,2)*XW(2,1)
        DDN(1,3)= Y12
        DDN(2,3)=-X12
C

        DO I=1,NODE
          DO J=1,NEQ
            DDN(J,I)=DDN(J,I)/DJAC
          ENDDO
        ENDDO
C
C
C ===== FOR THE CASE OF A 3-D TETRAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3)THEN
        DJAC=0.0D0
        DO KK=1,4
          IF(KK.EQ.1)THEN
            K1=2
            K2=3
            K3=4
          ELSEIF(KK.EQ.2)THEN
            K1=1
            K2=3
            K3=4
          ELSEIF(KK.EQ.3)THEN
            K1=1
            K2=2
            K3=4
          ELSE
            K1=1
            K2=2
            K3=3
          ENDIF
C
          A(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(2,K2)*XW(3,K3)+
     >          XW(2,K1)*XW(3,K2)*XW(1,K3)+XW(3,K1)*XW(1,K2)*XW(2,K3)-
     >          XW(1,K3)*XW(2,K2)*XW(3,K1)-XW(2,K3)*XW(3,K2)*XW(1,K1)-
     >          XW(3,K3)*XW(1,K2)*XW(2,K1))
C
          B(KK)=(-1.0D0)**KK*(XW(2,K1)*XW(3,K2)+
     >          XW(2,K2)*XW(3,K3)+XW(2,K3)*XW(3,K1)-
     >          XW(2,K3)*XW(3,K2)-XW(2,K2)*XW(3,K1)-
     >          XW(2,K1)*XW(3,K3))
C
          C(KK)=(-1.0D0)**(KK+1)*(XW(1,K1)*XW(3,K2)+
     >          XW(1,K2)*XW(3,K3)+XW(1,K3)*XW(3,K1)-
     >          XW(1,K3)*XW(3,K2)-XW(1,K2)*XW(3,K1)-
     >          XW(1,K1)*XW(3,K3))
C
          D(KK)=(-1.0D0)**KK*(XW(1,K1)*XW(2,K2)+
     >          XW(1,K2)*XW(2,K3)+XW(1,K3)*XW(2,K1)-
     >          XW(1,K3)*XW(2,K2)-XW(1,K2)*XW(2,K1)-
     >          XW(1,K1)*XW(2,K3))
C
          DJAC=DJAC+A(KK)
C          DN(KK)=A(KK)+B(KK)*XTEMP(1)+C(KK)*XTEMP(2)+D(KK)*XTEMP(3)
          DDN(1,KK)=B(KK)
          DDN(2,KK)=C(KK)
          DDN(3,KK)=D(KK)
        ENDDO
C MWF REPLACE WITH ABS IN CASE HAVE NEGATIVE JACOBIANS?
C        DJAC=DABS(DJAC)
        DO KK=1,4
          DO J=1,NEQ
            DDN(J,KK)=DDN(J,KK)/DJAC
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Start working on analytical tracking in separate routine, then merge
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PT123A
     I     (IDVE,IDVT,IVERBOSE,MAXEQ,
     I      NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O,
     I      ATOL,RTOL,SF, DN_SAFE,
     I      XG,IE,NLRL,LRL,IB,NORMALS,XFBAR,
     I      VTL2G,VT1E,VT2E,T1,T2,
     I      X_I,T_I,
     M      TPT,MPT,IDPT,
     O      XPT)
     
C 
C 05/15/2010 (MWF)
C ======================================================================
C < PURPOSE > 
C   IMPLEMENT PARTICLE TRACKING USING SEMI-ANALYTICAL TRACKING ON AN ELEMENT-BY-
C   ELEMENT BASIS IN THE DESIGNATED UNSTRUTURED MESH
C     ***NOTE***
C     EXTERNAL CODE MUST CHANGE CONNECTIVITY, LOCAL2GLOBAL, AND FLAG ARRAYS
C       TO MAKE SURE THEY ARE CONSISTENT WITH BASE ZERO 
C     THIS ROUTINE INTERNALLY SCALES THE VELOCITY BY +/- 1 BASED ON DIR
C
C
C
C
C
C TODO 
C      
C ======================================================================
C
      IMPLICIT NONE
C --HOW MUCH TO PRINT OUT 
      INTEGER IVERBOSE
Cf2py integer optional, intent (in) :: IVERBOSE = 0
C --MESH REPRESENTATION--
C THESE COULD BE PASSED IN TO MAKE DIMENSIONING SIMPLER
      INTEGER MAXEQ,MAXND
      PARAMETER(MAXND=8)
Cf2py integer optional, intent (in) :: MAXEQ = 3
C nNodes_global,nElements_global,nNodes_element,nSpace,nPointsToTrack,
C direction, output stream id
      INTEGER NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O
Cf2py integer required, intent (in) :: NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O
C NODE COORDS (nodeArray)
      DOUBLE PRECISION XG(MAXEQ*NNP)
Cf2py  double precision, intent (c), intent (in) :: XG(MAXEQ*NNP)
C ELEMENT NODE LOOKUP (elementNodesArray)
      INTEGER IE(NNDE*NEL)
Cf2py integer, intent (in) :: IE(NNDE*NEL)
C NODE - ELEMENTS IN NODE STAR LOOKUP (nodeElementOffsets,nodeElementsArray)
      INTEGER LRL(*),NLRL(*)
Cf2py integer, intent (in)  :: LRL(*),NLRL(*)
C FLAG ARRAY TO MARK NODES THAT ARE ON EXTERIOR BOUNDARY
      INTEGER IB(NNP)
Cf2py integer, intent (in)  :: IB(NNP)
C ELEMENT BOUNDARY UNIT OUTER NORMALS 
      DOUBLE PRECISION NORMALS(NEQ*NNDE*NEL)
Cf2py double precision, intent (in) :: NORMALS(NEQ*NNDE*NEL)
C ELEMENT BOUNDARY BARYCENTERS
      DOUBLE PRECISION XFBAR(MAXEQ*NNDE*NEL)
Cf2py double precision, intent (in) :: XFBAR(MAXEQ*NNDE*NEL)

C --TRACKING PARAMETERS--
      DOUBLE PRECISION ATOL,RTOL,SF,DN_SAFE
Cf2py double precision, intent (in) :: ATOL,RTOL,SF,DN_SAFE
C --TRACKING VELOCITY REPRESENTATION--
C TYPE OF LOCAL VELOCITY REPRESENTATION
C 1 -- 2 LOCAL SPACE IS C^0, P^1 
C 1 -- ASSUMED GLOBALLY CONTINUOUS (NODAL REPRESENTATION)
C 2 -- MAY BE DISCONTINUOUS (ELEMENT-BASED REPRESENTATION)
C      THIS VERSION DOES'T REALLY DISTINGUSIH SINCE 
C      VTL2G SHOULD ACCOUNT FOR THIS
C 3 -- RT0 WITH LOCAL BASIS \vec N_i = \frac{1}{d|\Omega_e|}(\vec x - x_{n,i}),
C 4 -- RT0 WITH LOCAL BASIS \vec N_i = \vec e_i i=0,...,d-1 and \vec N_d = \vec x
C SEE EL_VEL_PREP FOR DETAILS
      INTEGER IDVE
Cf2py integer, intent (in) :: IDVE = 2
C TYPE OF TIME DEPENDENCE IN VELOCITY FIELD
C 0 STEADY STATE
C 1 STRICTLY LINEAR VELOCITY DEPENDENCE (SEE BELOW)
      INTEGER IDVT
Cf2py integer, intent (in) :: IDVT = 0
C SEPARATE LOCAL TO GLOBAL MAP, 
C      INTEGER VTL2G(NEQ*NNDE*NEL)
CCf2py integer, intent (in) :: VTL2G(NEQ*NNDE*NEL)
      INTEGER VTL2G(*)
Cf2py integer, intent (in) :: VTL2G(*)

C DEGREES OF FREEDOM (cvelocity_dof)
      DOUBLE PRECISION VT1E(*),VT2E(*)
Cf2py double precision, intent (in) :: VT1E(*),VT2E(*)
C TIME LEVELS FOR VELOCITIES VT1E,VT2E
      DOUBLE PRECISION T1,T2
Cf2py double precision, intent (in) :: T1,T2
C

C -- INPUT TRACKING POINT INFORMATION
C POINTS TO TRACK (x_in) AND THEIR INITIAL TIMES (x_depart_times)
      DOUBLE PRECISION X_I(MAXEQ*NPT),T_I(NPT)
Cf2py double precision, intent (in) :: X_I(MAXEQ*NPT),T_I(NPT)
C
C BASE ONE, IDPT INPUT IS
C IF IDPT(K) >= 1 MEANS X(K) IS MESH NODE IDPT(K), TREAT DIFFERENTLY
C IF IDPT(K) == 0 THEN X(K) IS AN INTERIOR POINT
C IF IDPT(K) <  0 THEN DO NOT TRACK THE POINT
C ON OUTPUT SHOULD BE (BASE ONE)
C                   0  INTERIOR
C                  -1  EXITED DOMAIN SOMEWHERE 
C                  -2  DID NOT TRACK
      INTEGER IDPT(NPT)
Cf2py integer, intent (inplace) :: IDPT(NPT)
C -- IN/OUT TRACKING INFORMATION --
C TIMES TO TRACK POINTS TO (IN) FINAL TIME REACHED (OUT) 
      DOUBLE PRECISION TPT(NPT)
Cf2py double precision, intent (inplace) :: TPT(NPT)
C ELEMENT LOCATIONS FOR TRACKED POINTS
      INTEGER MPT(NPT)
Cf2py integer, intent (inplace) :: MPT(NPT)
C -- TRACKING OUTPUT
C LOCATION OF TRACKED POINTS AT END OF TRACKING CALL
      DOUBLE PRECISION XPT(MAXEQ*NPT)
Cf2py double precision, intent (inplace) :: XPT(MAXEQ*NPT)
C
C -- LOCAL VARIABLES --
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
      DOUBLE PRECISION DN_S(MAXND),DN(MAXND),DDN(MAXEQ,MAXND)
C
      DOUBLE PRECISION XOUT5(MAXEQ),XOUT4(MAXEQ),XERR(MAXEQ)
      DOUBLE PRECISION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ)
      DOUBLE PRECISION AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ)
C 
      DOUBLE PRECISION XS(MAXEQ),XTEMP(MAXEQ),VTEMP(MAXEQ)
      DOUBLE PRECISION DEQ,TS,DTS,T,SDT,DT0
      DOUBLE PRECISION DIR
C     
      INTEGER I,K,IPT,NODE,NP,M
      INTEGER IDSDT,I1,I2,I3,M2,M3,N1,N2,N3,J1,J2,J3
C MAKE THIS INPUT OR JUST GRAB FROM ATOL?
      DOUBLE PRECISION ZEROTOL
C
      ZEROTOL = ATOL
C
C
C ===== INITIALIZATION
C
      IF (MAXEQ.NE.3) THEN
         WRITE(LU_O,*)'pt123 currently requires MAXEQ=3, but MAXEQ= ',
     &        MAXEQ
         RETURN
      ENDIF
      NODE = NNDE
      DEQ  = 1.0E0/DBLE(NEQ)
      DIR = 1.D0
      IF (IBF.EQ.-1) THEN
         DIR = -1.D0
      ENDIF
C MWF COPY X_I INTO XPT
C     HAVE TO TREAT AS FLAT ARRAY, SINCE MAY BE DIFFERENT SHAPES IN PYTHON (q, versus ebqe)
C     ASSUME MPT,TPT SET BY CALLING CODE 
CMWF DEBUG
C      WRITE(LU_O,*)'ENTERING PT123A'
C      WRITE(LU_O,*)'NNP= ',NNP,' NEL= ',NEL,' NNDE= ',NNDE,' NEQ= ',NEQ,
C     &     ' NPT= ',NPT,' IBF= ',IBF,' LU_O= ',LU_O, ' ATOL= ',ATOL,
C     &     ' RTOL= ',RTOL, ' SF= ',SF, 'DN_SAFE= ',DN_SAFE

      DO IPT=1,NPT
CMWF DEBUG
C         WRITE(LU_O,*)'T_I(',IPT,')= ',T_I(IPT),
C     &        ' TPT(',IPT,')= ',TPT(IPT)
         DO I=1,NEQ
            XS(I)=X_I(I + MAXEQ*(IPT-1))
            XPT(I + MAXEQ*(IPT-1))=XS(I)
CMWF DEBUG
C            WRITE(LU_O,*)'X_I(',I,',',IPT,')= ',X_I(I+MAXEQ*(IPT-1)),
C     &           ' XPT(',I,',',IPT,')= ',XPT(I + MAXEQ*(IPT-1))

        ENDDO
      ENDDO
C
C =================== START PT USING ADAPTIVE RK ====================
C
C STEP 1.  LOCATE THE VERY FIRST TIME INTERVAL FOR PT
C
C
C USE T_START AS THE REFERENCE TIME FOR PT
C
C HAVE TO SET RRTIME FOR EACH POINT NOW
C
C 1.1 FOR THE CASE OF USING A STEADY VELOCITY FIELD
C
C
C STEP 2. CONDUCT PT WITHIN THE 1-ST TIME INTERVAL
C      
      DO 600 IPT=1,NPT
CMWF NOTE HAVE TO ACCOUNT FOR BASE ZERO INDEXING IN PROTEUS
C    AND MODIFIED CONVENTION FOR IDPT
        NP=IDPT(IPT)
C
C WHEN NP IS NEG, NO FURTHER TRACKING IS CONDUCTED
C
        IF(NP.LT.0)GOTO 600
C MWF NPATH NOT STORED FOR NOW
C        KPATH=NPATH(IPT)
        TS=T_I(IPT)
        DTS=TPT(IPT)-TS
        T=TS
        SDT=DTS
        DT0=SDT
        DT0=DMIN1(DT0,T2-T)
        I1=-1
        I2=-1
        I3=-1
C MWF DEBUG
C        WRITE(LU_O,*)'START OF STEP 2 IPT= ',IPT,' T1= ',T1,' T2= ',T2,
C     &       ' DTS= ',DTS,' NP= ',NP
CMWF NOW HAVE TO REFERENCE AS FLAT TO HANDLE DIFFERENT SHAPES FROM PYTHON
        DO I=1,NEQ
          XS(I)=XPT(I + MAXEQ*(IPT-1))
C MWF DEBUG
C          WRITE(LU_O,*)'XS(',I,')= ',XS(I)
        ENDDO
        
        IF(NP.EQ.0)GOTO 200
C
C ===== FOR THE CASE THAT MPT(IPT)=0, I.E., THE VERY FIRST
C     TRACKING: LOOP OVER ALL CONNECTED ELEMENTS
C
C MWF DEBUG
C        WRITE(LU_O,*)' ENTERING NODE TRACKING STEP IPT= ', IPT
        DO 150 I=NLRL(NP)+1,NLRL(NP+1)
          M=LRL(I)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)

C
C CONDUCT TRACKING WITHIN ELEMENT M
C
          CALL ELTRAK123A
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,
     I         IPT, T1,T2, DEQ, DN_SAFE,
     I         ZEROTOL,IDVT,
     I         XW,VT1W,VT2W, 
     I         NORMALS(NEQ*NNDE*(M-1)+1),XFBAR(MAXEQ*NNDE*(M-1)+1),
     M         T,DT0,SDT,XS,
     M         XTEMP,XOUT5,VTEMP,
     M         DN_S,DN,DDN,
     O         IDSDT,XPT,TPT,I1,I2,I3)
C        NPATH(IPT)=KPATH
C
          IF(IDSDT.EQ.-1)THEN
            MPT(IPT)=M
C
C 0 -- INTERIOR
            IDPT(IPT)=0
            GOTO 600
          ENDIF
          IF(IDSDT.EQ.1)GOTO 250
 150    CONTINUE
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
        IF (IB(NP).EQ.1) THEN
           IDPT(IPT) = -1
        ELSE
           IDPT(IPT) = -2
        ENDIF
        IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &     (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
          WRITE(LU_O,*)'WARNING (1) IN PT123!!!' 
          WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
          WRITE(LU_O,*)'IPT= ',IPT
          WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
          WRITE(LU_O,*)'NA = NLRL(NP)+1 = ',NLRL(NP)+1
          WRITE(LU_O,*)'NB = NLRL(NP+1) = ',NLRL(NP+1)
          WRITE(LU_O,*)'LRL(NA..NB) =',(LRL(I),
     &         I=NLRL(NP)+1,NLRL(NP+1))
          WRITE(LU_O,*)
        ENDIF
        GOTO 600
C
  200   CONTINUE
C
C ===== FOR THE CASE THAT PT STARTS WITHIN AN ELEMENT,
C       I.E., MPT(IPT)=M
C
        M=MPT(IPT)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
C MWF DEBUG
C        WRITE(LU_O,*)' B4 ELTRACK ELEMENT LOOP T= ',T,' TPT= ',TPT(IPT),
C     &       ' IPT= ',IPT,' M= ',M
C        WRITE(LU_O,*)' XS= ',(XS(I),I=1,NEQ)
        CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C$$$        DO J=1,NODE
C$$$           IEM=IE(J + NNDE*(M-1))
C$$$           DO K=1,NEQ
C$$$              XW(K,J)=XG(K + MAXEQ*(IEM-1))
C$$$CMWF SCALE LOCAL VELOCITY BY DIRECTION
C$$$              IVDOF = VTL2G(K + NEQ*(J-1) + NEQ*NNDE*(M-1))
C$$$              VT1W(K,J)=VT1E(IVDOF)*DIR
C$$$              VT2W(K,J)=VT2E(IVDOF)*DIR
C$$$CMWF DEBUG
C$$$C              WRITE(LU_O,*)' IEM= ',IEM,' IVDOF= ',IVDOF,
C$$$C     &             ' XW(',K,',',J,')= ',XW(K,J)
C$$$C              WRITE(LU_O,*)' VT1W(',K,',',J,')= ',VT1W(K,J)


C$$$           ENDDO
C$$$        ENDDO

C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C MWF REMOVE MAXPATH, KPATH ARGS
        CALL ELTRAK123A
     I      (MAXEQ,MAXND,NPT,NEQ,NODE,
     I       IPT, T1,T2, DEQ, DN_SAFE,
     I       ZEROTOL,IDVT,
     I       XW,VT1W,VT2W,
     I       NORMALS(NEQ*NNDE*(M-1)+1),XFBAR(MAXEQ*NNDE*(M-1)+1),
     M       T,DT0,SDT,XS,
     M       XTEMP,XOUT5,VTEMP,
     M       DN_S,DN,DDN,
     O       IDSDT,XPT,TPT,I1,I2,I3)
CMWF        NPATH(IPT)=KPATH
C MWF DEBUG
C        WRITE(LU_O,*)' AFTER ELTRACK ELEMENT LOOP T= ',T,
C     &       ' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT,' M= ',M
C        WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
C        DO K=1,NEQ
C           WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',XPT(K + IPT*MAXEQ)
C        ENDDO
C
        IF(IDSDT.EQ.-1)THEN
          MPT(IPT)=M
C 0 -- INTERIOR
          IDPT(IPT)=0
          GOTO 600
        ENDIF
C
C ===== CONTINUE PT WITHIN A NEW ELEMENT 
C  
  250   CONTINUE
C MWF TRY TO CATCH SOMETHING WRONG WITH ELTRAK123
        IF (I1.EQ.-1.OR.I2.EQ.-1.OR.I3.EQ.-1) THEN
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
           WRITE(LU_O,*)'WARNING (5) IN PT123!!!' 
           WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
           WRITE(LU_O,*)'IPT= ',IPT
           WRITE(LU_O,*)'M= ',M
           WRITE(LU_O,*)'I1= ',I1,' I2= ',I2,' I3= ',I3
           WRITE(LU_O,*)'NP = IDPT(IPT) = ',IDPT(IPT)
C -2 -- FAILED
           IDPT(IPT) = -2
           GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON A TRIANGULAR FACE
C COMPOSED OF NODES N1, N2, N3, WHERE N1=IE(I1,M), N2=IE(I2,M), 
C N3=IE(I3,M)
C
        IF(I3.NE.0)THEN
C MWF I1,I2,I3 DETERMINED LOCALLY SO CAN BE BASE 1
          N1=IE(I1 + NNDE*(M-1))
          N2=IE(I2 + NNDE*(M-1))
          N3=IE(I3 + NNDE*(M-1))
          DO J1=NLRL(N1)+1,NLRL(N1+1)          
            M=LRL(J1)
            DO J2=NLRL(N2)+1,NLRL(N2+1)          
              M2=LRL(J2)
              DO J3=NLRL(N3)+1,NLRL(N3+1)          
                M3=LRL(J3)
                IF(M.EQ.M2 .AND. M2.EQ.M3)THEN
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
C ALSO, WHEN USING IDVE=2, KNOW LOCAL VELOCITY IS ELEMENT BASED
C IN SOME CASES WE MAY HAVE TO USE VELOCITY SPACE L2G HERE
                   CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &                  NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &                  XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C$$$                   DO J=1,NODE
C$$$                      IEM=IE(J + NNDE*(M-1))
C$$$                      DO K=1,NEQ
C$$$                         XW(K,J)=XG(K + MAXEQ*(IEM-1))
C$$$CMWF SCALE LOCAL VELOCITY BY DIRECTION
C$$$                         IVDOF = VTL2G(K + NEQ*(J-1) + NEQ*NNDE*(M-1))
C$$$                         VT1W(K,J)=VT1E(IVDOF)*DIR
C$$$                         VT2W(K,J)=VT2E(IVDOF)*DIR
C$$$                      ENDDO
C$$$                   ENDDO
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C MWF GET RID OF MAXPT,KPATH ARGS
                  CALL ELTRAK123A
     I             (MAXEQ,MAXND,NPT,NEQ,NODE,
     I              IPT, T1,T2, DEQ,DN_SAFE,
     I              ZEROTOL,IDVT,
     I              XW,VT1W,VT2W,
     I              NORMALS(NEQ*NNDE*(M-1)+1),XFBAR(MAXEQ*NNDE*(M-1)+1),
     M              T,DT0,SDT,XS,
     M              XTEMP,XOUT5,VTEMP,
     M              DN_S,DN,DDN,
     O              IDSDT,XPT,TPT,I1,I2,I3)
C MWF                   NPATH(IPT)=KPATH
C
C MWF DEBUG
C                  WRITE(LU_O,*)' AFTER I3.NE.0 ELTRACK ELEMENT LOOP T= '
C     &                 ,T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= '
C     &                 ,IPT,' M= ',M
C                  WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
C                  DO K=1,NEQ
C                     WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
C     &                    XPT(K + IPT*MAXEQ)
C                  ENDDO
                  IF(IDSDT.EQ.-1)THEN
                    MPT(IPT)=M
C 0 -- INTERIOR
                    IDPT(IPT)=0
                    GOTO 600
                  ENDIF
                  IF(IDSDT.EQ.1)GOTO 250
                ENDIF
              ENDDO
            ENDDO
          ENDDO
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1)THEN
            IDPT(IPT) = -1
          ELSE
            IDPT(IPT) = -2
          ENDIF
          IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &       (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
            WRITE(LU_O,*)'WARNING (2) IN PT123!!!' 
            WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
            WRITE(LU_O,*)'IPT = ',IPT
            WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
C     WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
            WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
            WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
            WRITE(LU_O,*)'N1, N2, N3 = ',N1,N2,N3
            WRITE(LU_O,*)
          ENDIF
          GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON AN EDGE COMPOSED OF 
C NODES N1 AND N2, WHERE N1=IE(I1,M), N2=IE(I2,M)
C
        IF(I2.NE.0)THEN
          N1=IE(I1 + NNDE*(M-1))
          N2=IE(I2 + NNDE*(M-1))
          DO J1=NLRL(N1)+1,NLRL(N1+1)          
            M=LRL(J1)
            DO J2=NLRL(N2)+1,NLRL(N2+1)          
              M2=LRL(J2)
              IF(M.EQ.M2)THEN
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
                 CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &                NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &                XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C$$$                 DO J=1,NODE
C$$$                    IEM=IE(J + NNDE*(M-1))
C$$$                    DO K=1,NEQ
C$$$                       XW(K,J)=XG(K + MAXEQ*(IEM-1))
C$$$C MWF SCALE LOCAL VELOCITY BY DIRECTION
C$$$                       IVDOF = VTL2G(K + NEQ*(J-1) + NEQ*NNDE*(M-1))
C$$$                       VT1W(K,J)=VT1E(IVDOF)*DIR
C$$$                       VT2W(K,J)=VT2E(IVDOF)*DIR
C$$$                    ENDDO
C$$$                 ENDDO
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
C REMOVE MAXPT,MAXPATH,KPATH ARGS
                CALL ELTRAK123A
     I            (MAXEQ,MAXND,NPT,NEQ,NODE,
     I             IPT, T1,T2, DEQ, DN_SAFE,
     I             ZEROTOL,IDVT,
     I             XW,VT1W,VT2W,
     I             NORMALS(NEQ*NNDE*(M-1)+1),XFBAR(MAXEQ*NNDE*(M-1)+1),
     M             T,DT0,SDT,XS,
     M             XTEMP,XOUT5,VTEMP,
     M             DN_S,DN,DDN,
     O             IDSDT,XPT,TPT,I1,I2,I3)
CMWF SKIP NPATH FOR NOW
C                NPATH(IPT)=KPATH
C
C MWF DEBUG
C                WRITE(LU_O,*)' AFTER I2.NE.0 ELTRACK ELEMENT LOOP T= '
C     &               ,T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT
C     &               ,' M= ',M
C                WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
C                DO K=1,NEQ
C                   WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
C     &                  XPT(K + IPT*MAXEQ)
C                ENDDO
                IF(IDSDT.EQ.-1)THEN
                  MPT(IPT)=M
C 0 -- INTERIOR
                  IDPT(IPT)=0
                  GOTO 600
                ENDIF
                IF(IDSDT.EQ.1)GOTO 250
              ENDIF
            ENDDO
          ENDDO
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1)THEN
            IDPT(IPT) = -1
          ELSE
            IDPT(IPT) = -2
          ENDIF
          IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &       (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
            WRITE(LU_O,*)'WARNING (3) IN PT123!!!' 
            WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
            WRITE(LU_O,*)'IPT = ',IPT
            WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
CMWF          WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
            WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
            WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
            WRITE(LU_O,*)'N1, N2 = ',N1,N2
            WRITE(LU_O,*)
          ENDIF
          GOTO 600
        ENDIF
C
C FOR THE CASE THAT THE NEW TRACKING STARTS ON GLOBAL NODE N1, 
C WHERE N1=IE(I1,M)
C
C === LOOP OVER ALL CONNECTED ELEMENTS
C
        N1=IE(I1 + NNDE*(M-1))
        DO J1=NLRL(N1)+1,NLRL(N1+1)          
          M=LRL(J1)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
C ALSO, WHEN USING IDVE=2, KNOW LOCAL VELOCITY IS ELEMENT BASED
C IN SOME CASES WE MAY HAVE TO USE VELOCITY SPACE L2G HERE
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C$$$          DO J=1,NODE
C$$$             IEM=IE(J + NNDE*(M-1))
C$$$             DO K=1,NEQ
C$$$                XW(K,J)=XG(K + MAXEQ*(IEM-1))
C$$$C MWF SCALE LOCAL VELOCITY BY DIRECTION
C$$$                IVDOF = VTL2G(K + NEQ*(J-1) + NEQ*NNDE*(M-1))
C$$$                VT1W(K,J)=VT1E(IVDOF)*DIR
C$$$                VT2W(K,J)=VT2E(IVDOF)*DIR
C$$$             ENDDO
C$$$          ENDDO
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
CMWF REMOVE MAXPT, MAXPATH ARGS, KPATH
          CALL ELTRAK123A
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,
     I         IPT, T1,T2, DEQ, DN_SAFE,
     I         ZEROTOL,IDVT,
     I         XW,VT1W,VT2W,
     I         NORMALS(NEQ*NNDE*(M-1)+1),XFBAR(MAXEQ*NNDE*(M-1)+1),
     M         T,DT0,SDT,XS,
     M         XTEMP,XOUT5,VTEMP,
     M         DN_S,DN,DDN,
     O         IDSDT,XPT,TPT,I1,I2,I3)
CMWF          NPATH(IPT)=KPATH
C
C MWF DEBUG
C          WRITE(LU_O,*)' AFTER I1.NE.0 ELTRACK ELEMENT LOOP T= '
C     &         ,T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT
C     &         ,' M= ',M
C          WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
C          DO K=1,NEQ
C             WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',XPT(K + IPT*MAXEQ)
C          ENDDO
          IF(IDSDT.EQ.-1)THEN
            MPT(IPT)=M
C 0 -- INTERIOR
            IDPT(IPT)=0
            GOTO 600
          ENDIF
          IF(IDSDT.EQ.1)GOTO 250
        ENDDO  
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
        IF(IB(N1).EQ.1) THEN
          IDPT(IPT) = -1
        ELSE
          IDPT(IPT)=-2
        ENDIF
        IF((IVERBOSE.GE.2.AND.IDPT(IPT).EQ.-2).OR.
     &     (IVERBOSE.GE.3.AND.IDPT(IPT).EQ.-1)) THEN
          WRITE(LU_O,*)'WARNING (4) IN PT123!!!' 
          WRITE(LU_O,*)'THE IPT-TH PARTICLE CANNOT BE TRACKED:'
          WRITE(LU_O,*)'IPT = ',IPT
          WRITE(LU_O,*)'IDPT(IPT) = ',IDPT(IPT)
C        WRITE(LU_O,*)'KPATH = NPATH(IPT) = ',NPATH(IPT)
          WRITE(LU_O,*)'XS(1..NEQ) =',(XS(I),I=1,NEQ)
          WRITE(LU_O,*)'TPT(IPT) = ',TPT(IPT)
          WRITE(LU_O,*)'N1 = ',N1
          WRITE(LU_O,*)
        ENDIF
C     
  600 CONTINUE
        
      RETURN
      END
C
C 
C
C
C 
C
      SUBROUTINE ELTRAK123A
     I    (MAXEQ,MAXND,MAXPT,NEQ,NODE,
     I     IPT, T1,T2, DEQ, DN_SAFE,
     I     ZEROTOL,IVFLAG,
     I     XW,VT1W,VT2W,
     I     FNORMALS,XFBAR,
     M     T,DT0,SDT,XS,
     M     XTEMP,XOUT5,VTEMP,DN_S,DN,DDN,
     O     IDSDT,XPT,TPT,I1,I2,I3)
C 
C 02/08/2010 (MWF)
C ======================================================================
C < PURPOSE > 
C IMPLEMENT ANALYTICAL PARTICLE TRACKING WITHIN AN ELEMENT 
C TRACKING IN DIRECTION DIR
C IVFLAG = 1 RT0 VELOCITY, STRICTLY LINEAR VELOCITY DEPENDENCE
C
C IVFLAG = 0 (DEFAULT) STEADY-STATE RT0 VELOCITY ON A SIMPLEX
C 
C RT0 VELOCITY CAN BE REPRESENTED AS (CAPITAL LETTERS ARE VECTORS DIM=NEQ)
C
C V = V_0 + vx (X-X_0)
C
C WHERE X_0 IS THE STARTING POSITION WITH VELOCITY V_0=V(X_0)
C THE ANALYTICAL SOLUTION FOR THE POSITION IS
C 
C X(t) = X_0 + a(t)V_0
C
C  a(t)= (exp(vx * dt) - 1)/vx, if vx != 0
C      = dt,                     otherwise
C
C POINT X IS IN SIMPLEX M IFF (X_F-X).N_F <= 0 FOR ALL F IN FACES(M)
C WHERE 
C  X_F IS THE BARYCENTER OF FACE F, WITH UNIT OUTER NORMAL N_F
C 
C TO PERFORM TRACKING, WE CALCULATE TIMES TO INTERSECT BOUNDARIES OF M
C  
C (X_F - X(t)).N_F = 0
C (X_F - X_0 - a(t)V_0).N_F = 0
C (X_F - X_0).N_F  = a(t)V_0.N_F
C
C THEN IF V_0 != 0
C  dt = dx_F/v_F,  IF vx = 0
C     = ln(vx*dx_F/v_F + 1.0)/vx, IF vx != 0 and vx*dx_F/v_F + 1.0 > 0
C     = infty, OTHERWISE (NO INTERSECTION)
C WHERE
C  dx_F = (X_F-X_0).N_F,  v_F= V_0.N_F 
C
C IF V_0 = 0, THEN THERE IS NO INTERSECTION
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
CMWF NOW HAVE TO REFERENCE XPT AS FLAT ARRAY 
      DIMENSION XPT(MAXEQ*MAXPT),TPT(MAXPT)
      DIMENSION FNORMALS(NEQ,NODE),XFBAR(MAXEQ,NODE)
      DIMENSION XS(MAXEQ)
      DIMENSION XTEMP(MAXEQ),XOUT5(MAXEQ),VTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND),DDN(MAXEQ,MAXND)
      DIMENSION ICHECK(8)
      
C--- LOCAL VARIABLES ---
C TIME STEP TO TAKE, TIME TO INTERSECT BOUNDARY,
C   V.N_F, (X-X_0).N_F, vx
      DOUBLE PRECISION DT,DIR
      INTEGER IDEBUG
C
C =================== START TRACKING ====================
C
C TRACKING IN ELEMENT M
C
      ICOUNT=0
      IDSDT=1
      DT=DT0
      DIR = 1.D0
      IF (DT0.LT.0.D0) THEN 
         DIR = -1.D0
      ENDIF
      IDEBUG = 0
      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'ELTRAK123A CALLING STEP IVFLAG= ',IVFLAG,
     &        ' DT= ',DT,' XS= ',(XS(I),I=1,NEQ)
         DO J=1,NODE
            WRITE(6,*)' ',J,' XW= ',(XW(II,J),II=1,NEQ)
         ENDDO
      ENDIF

      IF (IVFLAG.EQ.1) THEN
         CALL STEPDTRT0V1
     I        (MAXEQ,MAXND,NEQ,NODE,
     I        DEQ, DN_SAFE,
     I        ZEROTOL,DIR,T1,T2,
     I        XW,VT1W,VT2W,
     I        FNORMALS,XFBAR,
     M        T,DT,XS,XOUT5,XTEMP,
     M        VTEMP,DN_S,DN,DDN)
      ELSE
C STEADY-STATE RT0
         CALL STEPDTSSRT0
     I        (MAXEQ,MAXND,NEQ,NODE,
     I        DEQ, DN_SAFE,
     I        ZEROTOL,DIR,T1,T2,
     I        XW,VT1W,VT2W,
     I        FNORMALS,XFBAR,
     M        T,DT,XS,XOUT5,XTEMP,
     M        VTEMP,DN_S,DN,DDN)
      ENDIF

C
C === EXAMINE THE COMPUTED ENDING LOCATION
C
      XSI=0.0D0
      CALL INTRP123
     I     (MAXEQ,MAXND,NEQ,NODE, XOUT5, XW,
     O     DN,DJAC)
      DL=DABS(DJAC)**(DEQ)
      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'XOUT5= ',(XOUT5(I),I=1,NEQ)
         WRITE(6,*)'DJAC= ',DJAC, 'DN= ',(DN(I),I=1,NODE)
      ENDIF
C
C === ADJUST DN AND XOUT5 WHEN THE END LOCATION IS VERY CLOSE TO 
C     ELEMENT BOUNDARY
C
      JCOUNT=0
      JOUTSIDE=0
      DD=0.0E0
      DO J=1,NODE
         ICHECK(J)=0
         IF(DN(J).LT.0.0E0)THEN
            IF(DABS(DN(J)).LE.DN_SAFE)THEN
               DD=DD+DABS(DN(J))
               DN(J)=0.0E0
               ICHECK(J)=1
               JCOUNT=JCOUNT+1
            ELSE
               JOUTSIDE=JOUTSIDE+1
            ENDIF 
         ELSEIF(DN(J).GT.1.0E0)THEN
            IF(DABS(DN(J)-1.0E0).LE.DN_SAFE)THEN
               DD=DD+DABS(DN(J)-1.0E0)
               DN(J)=1.0E0
               ICHECK(J)=1
               JCOUNT=JCOUNT+1
            ELSE
               JOUTSIDE=JOUTSIDE+1
            ENDIF
         ENDIF
      ENDDO
C
      IF(JCOUNT.NE.0 .AND. JOUTSIDE.EQ.0)THEN
         DO J=1,NODE
            IF(ICHECK(J).EQ.0)THEN
               DN(J)=DN(J)-DD
               IF(DN(J).LT.0.0E0)THEN
                  DD=-DN(J)
                  DN(J)=0.0E0
               ELSEIF(DN(J).GT.1.0E0)THEN
                  DD=DN(J)-1.0E0
                  DN(J)=1.0E0
               ENDIF
            ENDIF
         ENDDO
         DO I=1,NEQ
            XOUT5(I)=0.0E0
            DO J=1,NODE
               XOUT5(I)=XOUT5(I)+DN(J)*XW(I,J)
            ENDDO
         ENDDO
      ENDIF                 
C
C === EXAMINE THE VALUES OF DN
C
      IF (IDEBUG.GT.1) THEN
CMWF DEBUG
         WRITE(6,*)' ELTRAK EXAMINE DN, DN_S VALUES'
         WRITE(6,*)'DN= ',(DN(II),II=1,NODE)
         WRITE(6,*)'DN_S= ',(DN_S(II),II=1,NODE)
      ENDIF
      DO 150 I=1,NODE

C
C A. WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C    ===> COMPUTE XSI
C
         IF(DN(I).LT.0.0E0)THEN
            IF (IDEBUG.GT.0) THEN
C MWF DEBUG
               WRITE(6,*)' DN(',I,')= ',DN(I),' OUTSIDE ELEMENT'
            ENDIF
            D1=DN_S(I)
            D2=DN(I)
            D12=DN_S(I)-DN(I)
C
C A1. WHEN THERE IS NO TRACKIING THROUGH THIS ELEMENT
C
            IF(DABS(D1*DL).LT.ZEROTOL)THEN
               IF (IDEBUG.GT.0) THEN
C MWF DEBUG
                  WRITE(6,*)' ELTRAK CASE A1 D1= ',D1,' D2= ',D2,
     &                 ' DL= ',DL
               ENDIF
C
C IF THE PT IS LEAVING THE ELEMENT
C ==> SET IDSDT TO 0, AND IDENTIFY I1, I2, I3 TO
C     CONDUCT PT IN AN ADJACENT ELEMENT
C
               IF(DABS(D2*DL).GT.ZEROTOL)THEN
                  CALL DN_CHECK
     I                 (MAXND,NODE,NEQ, DL,ZEROTOL,
     I                 DN_S,
     O                 I1,I2,I3)
                  IDSDT=0
                  IF (IDEBUG.GT.0) THEN
C MWF DEBUG
                     WRITE(6,*)' ELTRAK DN_CHECK LEAVING IDSDT= ',IDSDT,
     &                    ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
     &                    ' D1= ',D1, ' D2= ',D2,' D12= ',D12
                  ENDIF
                  RETURN
               ENDIF
C
C A2. WHEN THERE IS A TRACKING THROUGH THIS ELEMENT
C
            ELSE
               IF (IDEBUG.GT.0) THEN
C MWF DEBUG
                  WRITE(6,*)'ELTRAK CASE A2 D1= ',D1,' D2= ',D2,
     &                 ' DL= ',DL
               ENDIF
C
C IF THE ENDING LOCATION CAN BE CONSIDERED ON THE ELEMENT BOUNDARY
C ==> NO NEED TO REDUCE TIMESTEP
C
               IF(DABS(D2*DL).LE.ZEROTOL) THEN
                  IF (IDEBUG.GT.0) THEN
C MWF DEBUG
                     WRITE(6,*)'D2*DL.LE.ZEROTOL GOTO 150'
                  ENDIF
                  GOTO 150
               ENDIF
C
C IF THE ENDING LOCATION IS TRULY OUTSIDE OF THE ELEMENT
C ==> COMPUTE AN INTERPOLATION FACTOR
C
C MWF THIS SHOULD NOT HAPPEN WITH ANALYTICAL TRACKING
               XSI=DMAX1(XSI,D12/D1)
               WRITE(6,*)'PROBLEM ELTRAK123A OVERSTEPPED ELEMENT XSI= ',
     &              XSI,' T= ',T,' DT= ',DT
               WRITE(6,*)' D12= ',D12,' D1= ',D1,' D2= ',D2,' DL= ',DL,
     &              'ZEROTOL= ',ZEROTOL,' DN_SAFE= ',DN_SAFE
               WRITE(6,*)'JCOUNT= ',JCOUNT,' JOUTSIDE= ',JOUTSIDE
               WRITE(6,*)'XS= ',(XS(II),II=1,NEQ)
               WRITE(6,*)'XOUT5= ',(XOUT5(II),II=1,NEQ)
               DO J=1,NODE
                  WRITE(6,*)' ',J,' XW= ',(XW(II,J),II=1,NEQ)
               ENDDO
               CALL EXIT(1)
            ENDIF
         ENDIF
 150  CONTINUE
C
C === IF XSI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C MWF THIS SHOULD NOT HAPPEN WITH ANALYTICAL TRACKING
C
      IF(XSI.GT.1.0E0)THEN
         WRITE(6,*)'ELTRAK123A XSI= ',XSI, 'SHOULD NOT BE > 1 HERE'
         CALL EXIT(1)
      ENDIF
C          DTT=DT/XSI
C          DT=DTT
C          GOTO 100
C        ENDIF
C
C B. WHEN THE ENDING LOCATION IS EITHER WITHIN THE ELEMENT OR 
C    ON THE ELEMENT BOUDNARY
C ==> UPDATE INFORMATION FOR THE SUCCESSIVE PT
C
      T=T+DT
CMWF        KPATH=KPATH+1
CMWF        TPT(KPATH,IPT)=T
      TPT(IPT)=T
      DO I=1,NEQ
         XS(I)=XOUT5(I)
CMWF          XPT(KPATH,I,IPT)=XOUT5(I)
CMWF NOW FLATTEN INDEXING
         XPT(I + MAXEQ*(IPT-1))=XOUT5(I)
      ENDDO   
      SDT=SDT-DT
      DT0=DT0-DT
C
C IF THE TRACKING TIME IS COMPLETELY CONSUMED
C ==> SET IDSDT TO -1
C
      IF(DABS(SDT).LE.1.0E-10)THEN
         IDSDT=-1
         SDT=0.0E0
      ENDIF
      IF(DABS(DT0).LE.1.0E-10)THEN
         IDSDT=-1
         DT0=0.0E0
      ENDIF
C
C IF THE ENDING LOCATION IS ON THE ELEMENT BOUNDARY
C ==> EXIT PT IN THIS ELEMENT
C
      DO I=1,NODE
         IF (IDEBUG.GT.0) THEN
CMWF DEBUG
            WRITE(6,*)'ELTRAK AFTER 150 DN(',I,')= ',DN(I),
     &           ' DABS(DN*DL) < ZEROTOL= ',DABS(DN(I)*DL).LT.ZEROTOL
         ENDIF
         IF(DABS(DN(I)*DL).LT.ZEROTOL)THEN
            CALL DN_CHECK
     I           (MAXND,NODE,NEQ, DL,ZEROTOL,
     I           DN,
     O           I1,I2,I3)
            IF (IDEBUG.GT.0) THEN
C MWF DEBUG
               WRITE(6,*)' ELTRAK DN_CHECK ON BNDY RETUR IDSDT= ',IDSDT, 
     &              ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
     &              ' DN_S= ',(DN_S(II),II=1,NODE)
            ENDIF
            RETURN
         ENDIF
      ENDDO
C
C UPDATE THE TRACKING TIMESTEP (I.E., DT) FOR THE SUCCESSIVE TRACKING
C WITHIN THIS SAME ELEMENT
C
CMWF ORIG
C        IF(IDSDT.EQ.-1)RETURN
C MWF DEBUG
      IF(IDSDT.EQ.-1) THEN
         IF (IDEBUG.GT.0) THEN
C MWF DEBUG
            WRITE(6,*)' ELTRAK DN_CHECK ON BNDY RETURN IDSDT= ',IDSDT, 
     &           ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
     &           ' DN_S= ',(DN_S(II),II=1,NODE)
         ENDIF
         RETURN
      ENDIF
C MWF END DEBUG
C MWF THIS SHOULD NOT HAPPEN WITH ANALYTICAL TRACKING?
      WRITE(6,*)'PROBLEM ELTRAK123A REACHED END WITHOUT CONCLUSION= ',
     &     ' T= ',T,' DT= ',DT
      WRITE(6,*)'XS= ',(XS(II),II=1,NEQ)
      WRITE(6,*)'XOUT5= ',(XOUT5(II),II=1,NEQ)
      CALL EXIT(1)
C     DT=DT0
C        DT=DMIN1(DT,DTT)
C        ICOUNT=0
C        GOTO 100
C
      
C 
C  999 CONTINUE
      RETURN
      END

      SUBROUTINE STEPDTSSRT0
     I    (MAXEQ,MAXND,NEQ,NODE,
     I     DEQ, DN_SAFE,
     I     ZEROTOL,DIR,TV1,TV2,
     I     XW,VT1W,VT2W,
     I     FNORMALS,XFBAR,
     M     T,DT,XS,XOUT,XTEMP,
     M     VTEMP,DN_S,DN,DDN)
      IMPLICIT NONE
C 
C ======================================================================
C < PURPOSE > 
C GIVEN INPUT POSITION XS, TARGET TIME STEP DT, VELOCITY REP IN V 
C COMPUTE NEW POSITION XOUT, AND TIME INTERVAL DT FOR X 
C IN CURRENT ELEMENT
C TRACKING IN DIRECTION DIR
C  STEADY-STATE RT0 VELOCITY ON A SIMPLEX
C 
C RT0 VELOCITY CAN BE REPRESENTED AS (CAPITAL LETTERS ARE VECTORS DIM=NEQ)
C
C V = V_0 + vx (X-X_0)
C
C WHERE X_0 IS THE STARTING POSITION WITH VELOCITY V_0=V(X_0)
C THE ANALYTICAL SOLUTION FOR THE POSITION IS
C 
C X(t) = X_0 + a(t)V_0
C
C  a(t)= (exp(vx * dt) - 1)/vx, if vx != 0
C      = dt,                     otherwise
C
C POINT X IS IN SIMPLEX M IFF (X_F-X).N_F <= 0 FOR ALL F IN FACES(M)
C WHERE 
C  X_F IS THE BARYCENTER OF FACE F, WITH UNIT OUTER NORMAL N_F
C 
C TO PERFORM TRACKING, WE CALCULATE TIMES TO INTERSECT BOUNDARIES OF M
C  
C (X_F - X(t)).N_F = 0
C (X_F - X_0 - a(t)V_0).N_F = 0
C (X_F - X_0).N_F  = a(t)V_0.N_F
C
C THEN IF V_0 != 0
C  dt = dx_F/v_F,  IF vx = 0
C     = ln(vx*dx_F/v_F + 1.0)/vx, IF vx != 0 and vx*dx_F/v_F + 1.0 > 0
C     = infty, OTHERWISE (NO INTERSECTION)
C WHERE
C  dx_F = (X_F-X_0).N_F,  v_F= V_0.N_F 
C
C IF V_0 = 0, THEN THERE IS NO INTERSECTION
C ======================================================================
C NUMBER OF EQUATIONS (SPACE DIM), NODES, AND ACTUAL NUMBER
      INTEGER MAXEQ,MAXND,NEQ,NODE
C 1/NEQ, TOLERANCE FOR NEAR BOUNDARY, ZERO VELOCITY
      DOUBLE PRECISION DEQ, DN_SAFE, ZEROTOL
C DIRECTION INTEGRATING IN TIME (1.D0 OR -1.D0)
      DOUBLE PRECISION DIR
C TIME LEVELS AT WHICH HAVE VELOCITY INFORMATION
      DOUBLE PRECISION TV1,TV2
C CURRENT ELEMENTS NODES, NODAL VELOCITY REPRESENTATION
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
C CURRENT ELEMENTS OUTER NORMALS, FACE BARYCENTERS
      DOUBLE PRECISION FNORMALS(NEQ,NODE),XFBAR(MAXEQ,NODE)
C STARTING TIME, TARGET TIME STEP, STARTING POSITION, 
      DOUBLE PRECISION T,DT,XS(MAXEQ),XOUT(MAXEQ)
C TEMPORARY VARIABLE FOR POSITION, VELOCITY, 
      DOUBLE PRECISION XTEMP(MAXEQ),VTEMP(MAXEQ)
C SHAPE FUNCTIONS AT STARTING POSITION, CURRENT POSITION, AND GRADIENTS
      DOUBLE PRECISION DN_S(MAXND),DN(MAXND),DDN(MAXEQ,MAXND)
C LOCAL VARIABLES
      INTEGER I,J,II,IDEBUG,NFACE
      DOUBLE PRECISION DTF,VF,DXF,VX,DXDVF,DALPHA,DJAC,ETA,TVEVAL
C TOLERANCE FOR NEAR BOUNDARYNESS
      DOUBLE PRECISION BNDTOL,DL
      LOGICAL INTERSECT
C
      IDEBUG = 0
C ASSUMES NODE=NFACE=NEQ+1
      NFACE=NODE
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTSSRT0 BEFORE FACE LOOP T= ',T,' DT= ',DT,
     &        ' DIR= ',DIR
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
         DO J=1,NODE
            WRITE(6,*)'XW(',J,')= ',(XW(I,J),I=1,NEQ)
            WRITE(6,*)'FNORMALS(',J,')= ',(FNORMALS(I,J),I=1,NEQ)
         ENDDO
      ENDIF
C 
C 
C

C EVALUATE VELOCITY AND DERIVATIVE AT (X_0,T_0) = (XS,T)
C SINCE THIS ASSUMES VELOCITY IS CONSTANT OVER TIME STEP
C USE VELOCITY AT T+DT/2 AS V(*,T0)
C COULD USE T OR T+DT TOO
      TVEVAL = T+0.5*DT
      IF(DABS(TV2-TV1).LE.1.0D-9) THEN
         ETA=0.D0
      ELSE
         ETA=(TVEVAL-TV1)/(TV2-TV1)
      ENDIF
      CALL INTRP123
     I     (MAXEQ,MAXND,NEQ,NODE,XS,XW,
     O     DN,DJAC)
C USE DL TO SCALE TOLERANCES FOR BOUNDARIES
      DL = DABS(DJAC)**(DEQ)
C SET BNDTOL TO BE DIMENSIONAL BASED ON DN_SAFE WHICH IS USED FOR ELTRAK123A TESTS FOR LOCATION
C BASED ON SHAPE FUNCTIONS 
      BNDTOL = 0.1D0*DL*DN_SAFE
C SAVE LOCATION OF STARTING POINT
      DO I=1,NODE
         DN_S(I)=DN(I)
      ENDDO
      CALL DERIV123
     I     (MAXEQ,MAXND,NEQ,NODE,XS,XW,
     O     DDN,DJAC)
      DO I=1,NEQ
         VTEMP(I) = 0.D0
      ENDDO
      VX = 0.D0
      DO J=1,NODE
C ASSUMES RT0 SO DV_XX = DV_YY = DV_ZZ AND DV_XY=DV_YX=0
         IF (IDEBUG.GT.0) THEN
C MWF DEBUG
            WRITE(6,*)' STEPDTSSRT0 VX= ',VX,' VTAVG(1,',J,')= ',
     &           0.5*(VT1W(1,J)+VT2W(1,J)), 'DDN(1,',J,')= ',DDN(1,J)
         ENDIF
C         VX = VX + (1.D0-ETA)*VT1W(1,J)*DDN(1,J)+ETA*VT2W(1,J)*DDN(1,J)
         DO I=1,NEQ
C USE AVERAGE OF DV_XX,DV_YY,DV_ZZ
            VX = VX + DEQ*((1.D0-ETA)*VT1W(I,J)*DDN(I,J)+
     &           ETA*VT2W(I,J)*DDN(I,J))
C            VTEMP(I) = VTEMP(I)+0.5*(VT1W(I,J)+VT2W(I,J))*DN(J)
           VTEMP(I) = VTEMP(I)+(1.D0-ETA)*VT1W(I,J)*DN(J) + 
     &           ETA*VT2W(I,J)*DN(J)
         ENDDO
      ENDDO
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTSSRT0 VX= ',VX,' V0= ',(VTEMP(II),II=1,NEQ)
      ENDIF
C
C LOOP THROUGH FACES AND COMPUTE EXIT TIMES 
C
      DO J=1,NFACE
         INTERSECT = .TRUE.
         VF = 0.D0
         DXF= 0.D0
C INITIALIZE DTF EVEN THOUGH INTERSECT LOGIC SHOULD MAKE SURE NOT TESTED UNLESS IT'S 
C CALCULATED
         DTF= -DIR*1.0D32
         DO I=1,NEQ
            VF = VF + VTEMP(I)*FNORMALS(I,J)
            DXF= DXF+ (XFBAR(I,J)-XS(I))*FNORMALS(I,J)
         ENDDO
         IF (IDEBUG.GT.0) THEN
CMWF DEBUG
            WRITE(6,*)'STEPDTSSRT0 FACE= ',J,' DIR= ',DIR,
     &           ' VF= ',VF,' DXF= ',DXF
         ENDIF
C TODO, SWITCH ZEROTOL -> DN_SAFE?
         IF (DABS(VF).LE.ZEROTOL.OR.DABS(DXF).LE.BNDTOL) THEN
            INTERSECT = .FALSE.
         ELSE
C FIRST ASSUME VX=0
            DXDVF = DXF/VF
            DTF   = DXDVF
            IF (DABS(VX).GT.ZEROTOL) THEN
               IF ((DXDVF*VX + 1.D0).LE.0.0) THEN
                  INTERSECT = .FALSE.
               ELSE
                  DTF = DLOG(DXDVF*VX + 1.D0)/VX
               ENDIF
            ENDIF
         ENDIF

C CHECK TO SEE IF TRAJECTORY INTERSECTS BOUNDARY AND IF THIS IS IN THE RIGHT DIRECTION
         IF (INTERSECT.AND.DIR*DTF.GT.ZEROTOL) THEN
            IF (IDEBUG.GT.0) THEN
CMWF DEBUG
               WRITE(6,*)'STEPDTSSRT0 INTERSECTED FACE= ',J,' DIR= ',
     &              DIR,' DTF= ',DTF,' VF= ',VF,' DXF= ',DXF,' DT= ',DT
            ENDIF
            IF (DABS(DTF).LT.DABS(DT)) THEN
               DT = DTF
            ENDIF
         ENDIF
C END FACE LOOP
      ENDDO

C NOW COMPUTE END LOCATION ASSUMING TOOK STEP OF SIZE DT
      DALPHA = DT
      IF (DABS(VX).GT.ZEROTOL) THEN
         DALPHA = (DEXP(VX*DT)-1.D0)/VX
      ENDIF
      DO I=1,NEQ
         XOUT(I) = XS(I) + DALPHA*VTEMP(I)
      ENDDO

      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'STEPDTSSRT0 AFTER FACE LOOP T= ',T,' DT= ',DT,
     &        ' DALPHA= ',DALPHA
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
         WRITE(6,*)'XOUT= ',(XOUT(I),I=1,NEQ)
      ENDIF
      
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C BEGIN ROUTINES FOR RT0 VELOCITY FIELD WITH
C JUST LINEAR TIME DEPENDENCY (NOT BILINEAR IN X-T)
C (CAPITAL LETTERS ARE VECTORS)
C  
C FOR SIMPLICITY ASSUME LOCAL ELEMENT REPRESENTATION IS 
C  V_e  = A + b X
C
C WITH FULL LINEAR INTERPOLATION IN TIME HAVE OVER [t_0,t_1]
C
C  V(X,t) = (A_0 + b_0X)(1 - dt/dt_1) + (A_1 + b_1 X)dt/dt_1
C  
C OR
C
C  V = V_0 + v_x (X-X_0) + V_t (t-t_0) + v_xt (X-X_0)(t-t_0)
C
C EVALUATE v_x  = \grad V_xx
C           V_t = \partial V / \partial t = V(t+dt)-V(t_0) 
C           v_xt= \partial v_x / \partial t
C
C   \grad V_ij = [b_0(1 - dt/dt_1) + b_f dt/dt_1] \delta ij
C   V_t = (A_1 - A_0)/dt_1 + (b_1-b_0)/dt_1 X
C   v_xt= (b_1 - b_0)/dt_1
C
C 
C 'SIMPLIFIED' LINEAR DEPENDENCY IS 
C   V = V_0 + v_x (X-X_0) + V_t (t-t_0)
C   I.E., b_0 = b_1
C
C UPDATE FORMULA IS
C  X(dt) = X_0 + alpha(dt) V_0 + beta(dt) V_t 
C
C  alpha(dt) = (exp(v_x*dt) - 1)/v_x, v_x != 0
C            = dt, otherwise
C  beta(dt)  = (exp(v_x*dt) - (1 + v_x*dt))/(v_x*v_x), v_x != 0
C            = dt*dt/2, otherwise
C 
C 
C AS IN STEADY-STATE CASE, COMPUTE INTERSECTION TIME WITH FACE f AS
C  (X(dt_i) - X_f).N_f = 0
C 
C ACTUALLY USE NEGATIVE OF STEADY-STATE EXAMPLE FOR SIMPLICITY
C 
C 
C ROUTINE FOR EVALUATING COEFFICIENTS FOR POSITION AS A FUNCTION OF VELOCITY
C AND TIME STEP  
      SUBROUTINE DXCOEFS1
     I     (ZEROTOL,DT,VX,
     O     ALPHA,DALPHA,
     O     BETA,DBETA)
      IMPLICIT NONE
      DOUBLE PRECISION ZEROTOL,DT,VX
      DOUBLE PRECISION ALPHA,BETA
      DOUBLE PRECISION DALPHA,DBETA
C
      DOUBLE PRECISION EVXDT
C
      IF(DABS(VX).GT.ZEROTOL) THEN
         EVXDT = DEXP(VX*DT)
         ALPHA = (EVXDT - 1.D0)/VX
         DALPHA= EVXDT
         BETA  = (EVXDT - (1.D0 + VX*DT))/(VX**2)
         DBETA = (EVXDT - 1.D0)/VX
      ELSE
         ALPHA = DT
         DALPHA= 1.D0
         BETA  = 0.5*DT**2
         DBETA = DT
      ENDIF
      RETURN
      END
C
C
C RESIDUAL ROUTINE FOR INTERSECTION WITH FACE f
C  (X(dt) - X_f).N_f = R
C  dx_0f + alpha(dt) v_0f + beta(dt) v_tf = R
C
C  dx_0f = (X_0-X_f).N_f
C  v0_f  = V_0.N_f, v_tf = V_t.N_f
C  alpha(dt) = (exp(v_x*dt) - 1)/v_x, v_x != 0
C            = dt_1, otherwise
C  beta(dt)  = (exp(v_x*dt) - (1 + v_x*dt))/(v_x*v_x), v_x != 0
C            = dt*dt/2, otherwise
C  alpha'(dt)= exp(v_x*dt), v_x != 0
C            = 1, otherwise
C  beta'(dt) = (exp(v_x*dt)-1)/v_x, v_x != 0
C            = dt, otherwise
C 
      SUBROUTINE RESDXF1
     I     (DT,
     O     R,DR,
     M     XPAR,IPAR)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION DT
Cf2py double precision, intent (in) :: DT
C R = X(dt).N_f - X_f.N_f 
      DOUBLE PRECISION R
Cf2py double precision, intent (out) :: R
C dR/d(dt)
      DOUBLE PRECISION DR
Cf2py double precision, intent (out) :: DR
C ARRAYS FOR PASSING PARAMETERS
      DOUBLE PRECISION XPAR(10)
Cf2py double precision, intent (inout) :: XPAR(10)
      INTEGER IPAR(10)
Cf2py integer, intent (inout) :: IPAR(10)
C
C XPAR(1) = dx_0f
C XPAR(2) = v_0f
C XPAR(3) = v_tf
C XPAR(4) = v_x
C XPAR(5) = epsilon
C LOCAL VARIABLES
      DOUBLE PRECISION DX0F,V0F,VTF,VX,EVXDT,ALPHA,DALPHA,BETA,DBETA
      DOUBLE PRECISION ZEROTOL
      
      DX0F= XPAR(1)
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)
      
      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)

      R = DX0F + ALPHA*V0F + BETA*VTF
      DR= DALPHA*V0F + DBETA*VTF

      RETURN
      END
C
C RESIDUAL ROUTINE FOR VELOCITY CRITICAL POINT
C  V(X(dt),dt).N_f = R
C  [V_0 + v_x (X(dt)-X_0) + V_t dt].N_f      = R 
C  v_f(dt) = v_0f + v_x dx_f(dt) + v_tf dt   = R 
C  dx_f(dt) = alpha(dt)v_0f + beta(dt) v_tf
C  v_0f = V_0 .N_f, v_tf = V_t . N_f
C  
C DERIVATIVE IS
C v_f'(dt) = v_x dx_f'(dt) + v_tf 
C          = vx v_f(dt) 
C
C
C 
      SUBROUTINE RESVF1
     I     (DT,
     O     R,DR,
     M     XPAR,IPAR)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION DT
Cf2py double precision, intent (in) :: DT
C R = V(X(dt),dt).N_f
      DOUBLE PRECISION R
Cf2py double precision, intent (out) :: R
C dR/d(dt)
      DOUBLE PRECISION DR
Cf2py double precision, intent (out) :: DR
C ARRAYS FOR PASSING PARAMETERS
      DOUBLE PRECISION XPAR(10)
Cf2py double precision, intent (inout) :: XPAR(10)
      INTEGER IPAR(10)
Cf2py integer, intent (inout) :: IPAR(10)
C
C XPAR(1) = dx_0f
C XPAR(2) = v_0f
C XPAR(3) = v_tf
C XPAR(4) = v_x
C XPAR(5) = epsilon
C LOCAL VARIABLES
      DOUBLE PRECISION VF,DXF,V0F,VTF,VX,ALPHA,DALPHA,BETA,DBETA
      DOUBLE PRECISION ZEROTOL
            
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)

      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)

      DXF = ALPHA*V0F + BETA*VTF
      VF  = V0F + VX*DXF + VTF*DT
      R   = VF
      DR = VX*VF + VTF

      RETURN
      END

      SUBROUTINE STEPDTRT0V1
     I    (MAXEQ,MAXND,NEQ,NODE,
     I     DEQ, DN_SAFE,
     I     ZEROTOL,DIR,TV1,TV2,
     I     XW,VT1W,VT2W,
     I     FNORMALS,XFBAR,
     M     T,DT,XS,XOUT,XTEMP,
     M     VTEMP,DN_S,DN,DDN)
      IMPLICIT NONE
C 
C ======================================================================
C < PURPOSE > 
C GIVEN INPUT POSITION XS, TARGET TIME STEP DT, VELOCITY REP IN V 
C COMPUTE NEW POSITION XOUT, AND TIME INTERVAL DT FOR X 
C IN CURRENT ELEMENT
C TRACKING IN DIRECTION DIR
C
C TODO 
C  CHANGE VELOCITY INTERPOLATIN TO BE RELATIVE TO INPUT VELOCITY TIME LEVELS
C
C
C  RT0 VELOCITY ON A SIMPLEX WITH STRICTLY LINEAR VELOCITY DEPENDENCE (NO BILINEAR TERMS)
C 
C RT0 VELOCITY CAN BE REPRESENTED AS (CAPITAL LETTERS ARE VECTORS DIM=NEQ)
C
C  V = V_0 + v_x (X-X_0) + V_t (t-t_0) + v_xt (X-X_0)(t-t_0)
C
C WHERE X_0 IS THE STARTING POSITION WITH VELOCITY V_0=V(X_0)
C THE ANALYTICAL SOLUTION FOR THE POSITION IS
C 
C X(dt) = X_0 + alpha(dt)V_0 + beta(dt)V_t
C
C SEE ABOVE FOR alpha,beta formulas
C
C POINT X IS IN SIMPLEX M IFF (X_F-X).N_F <= 0 FOR ALL F IN FACES(M)
C WHERE 
C  X_F IS THE BARYCENTER OF FACE F, WITH UNIT OUTER NORMAL N_F
C 
C THIS APPROACH IS VERY SIMILAR TO THE STEADY-STATE CASE, EXCEPT THAT THERE
C CAN BE (AT MOST 1) A VELOCITY REVERSAL OVER A TIME STEP
C  
C TO PERFORM TRACKING, WE 
C  
C 1. CALCULATE FULL TIME STEP INFORMATION AS IF ELEMENT VELOCITY IS GLOBAL
C    a. X_1 = X_1(dt_1),  dt_1 = (t_1-t_0) INPUT TIME STEP
C    b. V_0 = V(X_0,t_0), V_1 = V(X_1,t_0+dt_1)
C 
C THEN FOR EACH FACE
C 2. COMPUTE 
C     dx_f0 = (X_0 - X_f).N_f, dx_f1 = (X_1-X_f).N_f, v_0f = V_0.N_f, v_1f = V_1.N_f
C
C    IF |dx_f0|  <= epsilon, SKIP FACE (ON BOUNDARY OR OUTSIDE ALREADY)
C
C    OTHERWISE EVALUATE CASES BELOW BASED ON v_0f, v_1f, WHERE dir IS THE 
C    DIRECTION GOING IN TIME, sign[v] IS CALCULATED USING v_t IF |v| <= 1.0e-3 epsilon
C    TODO CHECK THIS MODIFICATION OF HEALY AND RUSSELL, BUT I'M USING
C    sign[v_t*dir] IF |v_0f| <= 1.0e-3 epsilon and sign[-v_t*dir] IF |v_1f| <= 1.0e-3 epsilon
C    THEY DO NOT INCLUDE dir IN THE TEST
C
C A. sign[v_f0] = sign[v_f1] = dir
C    if dx_f1  < 0, never intersects face
C    otherwise,
C      solve for timestep  to intersect face, dt_i, using dx_f(dt) = R, bracketed between
C       R(0) = dx_f0 < 0, R(dt_1) = dx_f1 > 0   
C    set dt_out = min(dt_out,dt_i) 
C
C B. sign[v_f0] = sign[v_f1] = -dir
C    never intersects face, skip
C
C C. sign[v_f0] != sign[v_f1], sign[v_f0] = dir
C    solve for dt_s in [0,dt_1], with v_f(dt_s) = 0, using v_f(dt) = R, bracketed by v_0f, v_01
C    compute dx_s(dt_s)
C    if dx_s < 0, never intersects face
C    otherwise,
C      solve for timestep to intersect face, dt_i, dx_f(dt_i)= 0, using dx_f(dt) = R, bracketed between
C       R(0) = dx_f0 < 0, R(dt_s) = dx_s > 0,
C     set dt_out = min(dt_out,dt_i)
C
C D. sign[v_f0] != sign[v_f1], sign[v_f0] = -dir 
C    if dx_f1 < 0, never intersects
C    otherwise
C     solve for dt_s in [0,dt_1], with v_f(dt_s) = 0, using v_f(dt) = R, bracketed by v_0f, v_01
C     compute dx_s(dt_s)
C   
C    if dx_s > 0, should not be possible
C    otherwise,
C      solve for timestep to intersect face, dt_i, dx_f(dt_i)= 0, using dx_f(dt) = R, bracketed between
C       R(0) = dx_s < 0, R(dt_1) = dx_1 > 0,
C     set dt_out = min(dt_out,dt_i)
C ======================================================================
C NUMBER OF EQUATIONS (SPACE DIM), NODES, AND ACTUAL NUMBER
      INTEGER MAXEQ,MAXND,NEQ,NODE
C 1/NEQ, TOLERANCE FOR NEAR BOUNDARY, ZERO VELOCITY
      DOUBLE PRECISION DEQ, DN_SAFE, ZEROTOL
C DIRECTION INTEGRATING IN TIME (1.D0 OR -1.D0)
      DOUBLE PRECISION DIR
C TIME LEVELS AT WHICH VELOCITY DOFS ARE GIVEN  
      DOUBLE PRECISION TV1,TV2
C CURRENT ELEMENTS NODES, NODAL VELOCITY REPRESENTATION
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
C CURRENT ELEMENTS OUTER NORMALS, FACE BARYCENTERS
      DOUBLE PRECISION FNORMALS(NEQ,NODE),XFBAR(MAXEQ,NODE)
C STARTING TIME, TARGET TIME STEP, STARTING POSITION, 
      DOUBLE PRECISION T,DT,XS(MAXEQ),XOUT(MAXEQ)
C TEMPORARY VARIABLE FOR POSITION, VELOCITY, 
      DOUBLE PRECISION XTEMP(MAXEQ),VTEMP(MAXEQ)
C SHAPE FUNCTIONS AT STARTING POSITION, CURRENT POSITION, AND GRADIENTS
      DOUBLE PRECISION DN_S(MAXND),DN(MAXND),DDN(MAXEQ,MAXND)
C 
      DOUBLE PRECISION RTSAFE
      EXTERNAL RESVF1,RESDXF1
C LOCAL VARIABLES
      DOUBLE PRECISION X1(3),V1(3),VT(3),VTEMP2(3),RPAR(10)
      INTEGER I,J,II,IDEBUG,NFACE,IPAR(10)
      DOUBLE PRECISION DTF,VTF,DXF0,DXF1,VX,DTFI,DTFS,DTOUT,DXFS,DJAC,DL
      DOUBLE PRECISION V0F,V1F,ALPHA,DALPHA,BETA,DBETA,SV0F,SV1F,DXFI
      DOUBLE PRECISION ETAT0,ETAT1,ETATVX
C RES TOLERANCES FOR NONLINEAR SOLVER WHEN COMPUTING INTERSECTION, AND V=0
      DOUBLE PRECISION DXNLTOL,VNLTOL
C TOLERANCE FOR INCREMENT IN NONLINEAR SOLVER
      DOUBLE PRECISION DXINCRTOL,VINCRTOL
C TOLERANCE FOR NEAR BOUNDARYNESS
      DOUBLE PRECISION BNDTOL
C FOR NEWTON SOLVES
      DOUBLE PRECISION DTL,DTH
      INTEGER MAXIT
      LOGICAL INTERSECT
C
      IDEBUG = 0
C MWF HACK      
C      IF(T.GT.1.294) THEN
C         IDEBUG = 1
C      ENDIF
      DTOUT = DT
C ASSUMES NODE=NFACE=NEQ+1
      NFACE=NODE
C FOR NEWTON
      MAXIT = 100
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V1 BEFORE FACE LOOP T= ',T,' DT= ',DT,
     &        ' DIR= ',DIR
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
         DO J=1,NODE
            WRITE(6,*)'XW(',J,')= ',(XW(I,J),I=1,NEQ)
            WRITE(6,*)'XFBAR(',J,')= ',(XFBAR(I,J),I=1,NEQ)
            WRITE(6,*)'FNORMALS(',J,')= ',(FNORMALS(I,J),I=1,NEQ)
         ENDDO
      ENDIF
C 
C 
C
C EVALUATE VELOCITY AND DERIVATIVE AT (X_0,T_0) = (XS,T)
C ALLOW FOR VELOCITY DOFS TO BE GIVEN AT SOME TIME OTHER THAN
C T AND T+DT
      IF(DABS(TV2-TV1).LE.1.0D-9) THEN
         ETAT0 = 0.D0
         ETAT1 = 0.D0
         ETATVX= 0.D0
      ELSE
         ETAT0 = (T-TV1)/(TV2-TV1)
         ETAT1 = (T+DT-TV1)/(TV2-TV1)
         ETATVX= (T+DT*0.5D0-TV1)/(TV2-TV1)
      ENDIF
      CALL INTRP123
     I     (MAXEQ,MAXND,NEQ,NODE,XS,XW,
     O     DN,DJAC)
C USE DL TO SCALE TOLERANCES FOR BOUNDARIES
      DL = DABS(DJAC)**(DEQ)
C TIE TOLERANCE FOR BOUNDARY TO DN_SAFE USE BY CALLING ROUTINE, TRY TO BE A LITTLE TIGHTER TO
C ACCOUNT FOR ROUNDOFF ETC
C SET BNDTOL TO BE DIMENSIONAL BASED ON DN_SAFE WHICH IS USED FOR ELTRAK123A TESTS FOR LOCATION
C BASED ON SHAPE FUNCTIONS 
      BNDTOL = 0.1D0*DL*DN_SAFE
C      IF (DL.LT.0.1D0) BNDTOL = DL*DN_SAFE
C NEEDS TO BE CONSISTENT WITH NEAR BOUNDARY TEST BECAUSE RESIDUAL IS DISTANCE TO BOUNDARY
      DXNLTOL = BNDTOL
      DXINCRTOL= 0.1D0*DXNLTOL 
C MAKE CONSISTENT WITH TESTS WHEN VELOCITY IS POSTITIVE OR NEGATIVE
      VNLTOL  = ZEROTOL
      VINCRTOL= 0.1D0*VNLTOL 

C SAVE LOCATION OF STARTING POINT
      DO I=1,NODE
         DN_S(I)=DN(I)
      ENDDO
      CALL DERIV123
     I     (MAXEQ,MAXND,NEQ,NODE,XS,XW,
     O     DDN,DJAC)
      DO I=1,NEQ
         VTEMP(I) = 0.D0
         VTEMP2(I)= 0.D0
      ENDDO
      VX = 0.D0
      DO J=1,NODE
C ASSUMES RT0 SO DV_XX = DV_YY AND DV_XY=DV_YX=0
         IF (IDEBUG.GT.0) THEN
C MWF DEBUG
            WRITE(6,*)' STEPDTRT0V1 VX= ',VX,' VTAVG(1,',J,')= ',
     &           0.5*(VT1W(1,J)+VT2W(1,J)), 'DDN(1,',J,')= ',DDN(1,J)
         ENDIF
C NOW USE AVERAGE OF DV_XX, DV_YY, DV_ZZ
C         VX = VX + VT1W(1,J)*DDN(1,J)
         DO I=1,NEQ
            VX = VX + DEQ*((1.D0-ETATVX)*VT1W(I,J)*DDN(I,J)+
     &           ETATVX*VT2W(I,J)*DDN(I,J))
            VTEMP(I)  = VTEMP(I) + (1.D0-ETAT0)*VT1W(I,J)*DN(J) + 
     &           ETAT0*VT2W(I,J)*DN(J)
            VTEMP2(I) = VTEMP2(I)+ (1.D0-ETAT1)*VT1W(I,J)*DN(J) + 
     &           ETAT1*VT2W(I,J)*DN(J)
         ENDDO
      ENDDO
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V1 VX= ',VX,' V0= ',(VTEMP(II),II=1,NEQ),
     &        ' V(DT)= ',(VTEMP2(II),II=1,NEQ)
      ENDIF
C
C TIME DERIVATIVE, (CONSTANT IN SPACE), EVALUATED AT X_0 OVER FULL TIME STEP
C 
      DO I=1,NEQ
         VT(I) = (VTEMP2(I)-VTEMP(I))/DT
      ENDDO
C 
C STEP 1. COMPUTE POSITION AND VELOCITY AT END OF THE FULL TIME STEP
C
      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
      DO I=1,NEQ
         X1(I) = XS(I) + ALPHA*VTEMP(I) + BETA*VT(I)
         V1(I) = VTEMP(I) + VX*(X1(I)-XS(I)) + VT(I)*DT
      ENDDO
C
C
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V1 X1= ',(X1(II),II=1,NEQ),
     &        ' V(X1,DT)= ',(V1(II),II=1,NEQ)
      ENDIF

C
C LOOP THROUGH FACES AND COMPUTE EXIT TIMES 
C
      DO J=1,NFACE
         INTERSECT = .TRUE.
C INITIALIZE DTF EVEN THOUGH INTERSECT LOGIC SHOULD MAKE SURE NOT TESTED UNLESS IT'S 
C CALCULATED
         DTF= -DIR*1.0D32
C STEP 2.  
         V0F = 0.D0
         V1F = 0.D0
         DXF0= 0.D0
         DXF1= 0.D0
C INCLUDE V_t.N_f HERE FOR CONVENIENCE
         VTF = 0.D0
         DO I=1,NEQ
C NOTE DXF IS OPPOSITE SIGN OF STEADY-STATE EXAMPLE
            V0F = V0F + VTEMP(I)*FNORMALS(I,J)
            DXF0= DXF0+ (XS(I)-XFBAR(I,J))*FNORMALS(I,J)
            V1F = V1F + V1(I)*FNORMALS(I,J)
            DXF1= DXF1+ (X1(I)-XFBAR(I,J))*FNORMALS(I,J)
            VTF = VTF + VT(I)*FNORMALS(I,J)
         ENDDO
C SIGNS 
         IF (DABS(V0F).LE.1.D-3*ZEROTOL) THEN
C TODO HAVE TO CHECK THIS REASONING
            SV0F = DSIGN(1.D0,VTF*DIR)
         ELSE
            SV0F = DSIGN(1.D0,V0F)
         ENDIF
         IF (DABS(V1F).LE.1.D-3*ZEROTOL) THEN
C CHECK WHY THIS IS RIGHT ACCORDING TO RUSSELL AND HEALY
            SV1F = DSIGN(1.D0,-VTF*DIR)
         ELSE
            SV1F = DSIGN(1.D0,V1F)
         ENDIF

         IF (IDEBUG.GT.0) THEN
CMWF DEBUG
            WRITE(6,*)'STEPDTRT0V1 FACE= ',J,' DIR= ',DIR,
     &           'VTF= ',VTF,' V0F= ',V0F,' DXF0= ',DXF0,' SV0F= ',SV0F,
     &           ' V1F= ',V1F,' DXF1= ',DXF1,' SV1F= ',SV1F
         ENDIF
C TODO, SWITCH ZEROTOL -> BNDTOL?
C 
         IF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.GT.0.D0) THEN
C     CASE A
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE A DXF1.LT.-BNDTOL= ',DXF1.LT.-BNDTOL
            ENDIF
            IF (DXF1.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
               RPAR(1) = DXF0
               RPAR(2) = V0F
               RPAR(3) = VTF
               RPAR(4) = VX
               RPAR(5) = ZEROTOL
               DTL = 0.D0
               DTH = DT
               DTFI = RTSAFE(RESDXF1,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               DTF  = DTFI
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
                  CALL DXCOEFS1(ZEROTOL,DTFI,VX,ALPHA,DALPHA,BETA,DBETA)
                  DXFI = DXF0 + ALPHA*V0F + BETA*VTF
                  WRITE(6,*)'STEPDTRT0V1 CASE A DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFI= ',DTFI,' DTF= ',DTF,' DXFI= ',DXFI
                  WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &                 ' DL= ',DL,'IPAR(1)= ',IPAR(1),
     &                 ' IPAR(2)= ',IPAR(2)
               ENDIF
            ENDIF
         ELSEIF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.LT.0.D0) THEN
C     CASE B
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE B NO INTERSECTION'
            ENDIF
            INTERSECT = .FALSE.
         ELSEIF (SV0F*SV1F.LT.0.D0.AND.SV0F*DIR.GT.0.D0) THEN
C     CASE C
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE C DXF1.LT.-BNDTOL= ',DXF1.LT.-BNDTOL
            ENDIF
            RPAR(1) = DXF0
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX
            RPAR(5) = ZEROTOL
            DTL = 0.D0
            DTH = DT
            DTFS= RTSAFE(RESVF1,DTL,DTH,VINCRTOL,VNLTOL,MAXIT,RPAR,IPAR)
C     
            CALL DXCOEFS1(ZEROTOL,DTFS,VX,ALPHA,DALPHA,BETA,DBETA)
            DXFS = DXF0 + ALPHA*V0F + BETA*VTF 
            IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
               WRITE(6,*)'STEPDTRT0V1 CASE C DTL= ',DTL,' DTH= ',
     &              DTH,' DTFS= ',DTFS,' DXFS= ',DXFS
            ENDIF
            IF (DXFS.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
               DTL = 0.D0
               DTH = DTFS
               DTFI= RTSAFE(RESDXF1,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               DTF = DTFI
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
                  WRITE(6,*)'STEPDTRT0V1 CASE C DTL= ',DTL,
     &                 ' DTH= ',DTH,' DTFI= ',DTFI,
     &                 ' DTF= ',DTF
               ENDIF
            ENDIF
         ELSEIF (SV0F*SV1F.LT.0.D0.AND.SV0F*DIR.LT.0.D0) THEN
C     CASE D
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE D DXF1.LT.-BNDTOL= ',DXF1.LT.-BNDTOL
            ENDIF
            IF (DXF1.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
               RPAR(1) = DXF0
               RPAR(2) = V0F
               RPAR(3) = VTF
               RPAR(4) = VX
               RPAR(5) = ZEROTOL
               DTL = 0.D0
               DTH = DT
               DTFS= RTSAFE(RESVF1,DTL,DTH,VINCRTOL,VNLTOL,
     &              MAXIT,RPAR,IPAR)
C     
               CALL DXCOEFS1(ZEROTOL,DTFS,VX,ALPHA,DALPHA,BETA,DBETA)
               DXFS = DXF0 + ALPHA*V0F + BETA*VTF 
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
                  WRITE(6,*)'STEPDTRT0V1 CASE D DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFS= ',DTFS,' DXFS= ',DXFS
               ENDIF
               IF (DXFS.GT.BNDTOL) THEN
                  INTERSECT = .FALSE.
                  WRITE(6,*) 'PROBLEM CASE D DXFS = ',DXFS,' WRONG'
                  CALL EXIT(1)
               ELSE
                  DTL = DTFS
                  DTH = DT
                  DTFI= RTSAFE(RESDXF1,DTL,DTH,DXINCRTOL,DXNLTOL,
     &                 MAXIT,RPAR,IPAR)
                  DTF = DTFI
                  IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
                     WRITE(6,*)'STEPDTRT0V1 CASE D DTL= ',DTL,
     &                    ' DTH= ',DTH,' DTFI= ',DTFI,
     &                    ' DTF= ',DTF
                  ENDIF
               ENDIF
            ENDIF
C     END ALL CASES?
         ENDIF
C CHECK TO SEE IF TRAJECTORY INTERSECTS BOUNDARY WITH NONZERO STEP 
C AND IF THIS IS IN THE RIGHT DIRECTION
         IF (INTERSECT.AND.DIR*DTF.GT.ZEROTOL.AND.DABS(DTF).GT.ZEROTOL)
     &        THEN
            IF (IDEBUG.GT.0) THEN
CMWF DEBUG
               WRITE(6,*)'STEPDTRT0V1 INTERSECTED FACE= ',J,' DIR= ',
     &              DIR,' DTF= ',DTF,' V0F= ',V0F,' DXF0= ',DXF0,
     &              ' DTOUT= ',DTOUT
            ENDIF
            IF (DABS(DTF).LT.DABS(DTOUT)) THEN
               DTOUT = DTF
            ENDIF
         ENDIF
C END FACE LOOP
      ENDDO
C
      DT = DTOUT
C NOW COMPUTE END LOCATION ASSUMING TOOK STEP OF SIZE DT
      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
      
      DO I=1,NEQ
         XOUT(I) = XS(I) + ALPHA*VTEMP(I) + BETA*VT(I)
      ENDDO

      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'STEPDTRT0V1 AFTER FACE LOOP T= ',T,' DT= ',DT,
     &        ' ALPHA= ',ALPHA
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
         WRITE(6,*)'XOUT= ',(XOUT(I),I=1,NEQ)
      ENDIF

      RETURN
      END
C
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C NEWTON WITH BISECTION ROUTINES FOR LIU AND HEALY, RUSSELL ALGORITHMS
C TAKEN FROM NUMERICAL RECIPES IN FORTRAN
C
C Using a combination of Newton-Raphson method and bisection, 
C find the root of a function bracketed between x1 and x2.
C The root rtsafe will be refined until its accuracy is known within 
C +/- xacc. funcd  is a user-supplied subroutine that returns both the function value and 
C the first derivative  of the function at the point x. 
C add facc tolerance for residual and set ipar(1) to flag on convergence as well
      FUNCTION rtsafe(funcd,x1,x2,xacc,facc,maxit,xpar,ipar) 
      IMPLICIT NONE
      EXTERNAL funcd
C     MAX ITERATIONS
      INTEGER maxit
Cf2py intent (in) maxit :: maxit = 100
C ROOT      
      DOUBLE PRECISION rtsafe
C BRACKET INTERVAL, ERROR
      DOUBLE PRECISION x1,x2,xacc,facc
Cf2py intent (inout) x1,x2,xacc 
      DOUBLE PRECISION xpar(10)
      INTEGER ipar(10)
C local variables
      INTEGER j 
      DOUBLE PRECISION df,dx,dxold,f,fh,fl,temp,xh,xl
C -1 did not converge
C  0 converge on |dx| only
C  1 converge on |dx| and |facc| 
C  2 converge on |facc| only
      ipar(1) = -1
C number of iterations
      ipar(2) = 0
      call funcd(x1,fl,df,xpar,ipar)
      call funcd(x2,fh,df,xpar,ipar)
       
      if((fl.gt.facc.and.fh.gt.facc).or.
     &     (fl.lt.(-facc).and.fh.lt.(-facc))) then
         write(6,*) 'root must  be bracketed in rtsafe'
         call exit(1)
      endif

      if(dabs(fl).le.facc) then
         rtsafe = x1
         ipar(1) = 2
         return
      else if(dabs(fh).le.facc) then
         rtsafe = x2
         ipar(1)= 2
         return
      else if(fl.lt.0.D0) then
         xl=x1
         xh=x2
      else
         xh=x1
         xl=x2
      endif
C initialize guess for root
      rtsafe= 0.5*(x1+x2)
C stepsize before last
      dxold = dabs(x2-x1)
C and the last step
      dx=dxold

      call funcd(rtsafe,f,df,xpar,ipar)

      do j=1,maxit
         ipar(2) = j
         if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0.D0
     &        .or. dabs(2.D0*f).gt.dabs(dxold*df) ) then
            dxold=dx
            dx=0.5*(xh-xl)
            rtsafe=xl+dx
C change in root is negligible
            if(xl.eq.rtsafe) then
               ipar(1) = 0
C shouldn't be true
               if (dabs(fl).le.facc) ipar(1)=1
               return
            endif
         else
C Newton step acceptable
            dxold=dx
            dx=f/df
            temp=rtsafe
            rtsafe=rtsafe-dx
            if(temp.eq.rtsafe) then
               ipar(1) = 0
               if (dabs(f).le.facc) ipar(1)=1
               return
            endif
         endif
C 
         call funcd(rtsafe,f,df,xpar,ipar)
C convergence test, move after feval to use both tests
         if (dabs(dx).lt.xacc.and.dabs(f).le.facc) then
            ipar(1) = 1
            return
         endif
C maintain root bracketing
         if(f.lt.0.D0) then
            xl=rtsafe
         else
            xh=rtsafe
         endif
      enddo
      write(6,*)'rtsafe exceeded maximum its= ',maxit
      call exit(1)
      return
      end
      
