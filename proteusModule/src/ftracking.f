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
C CL      -- elementDiametersArray
C IB      -- nodeOnBoundaryArray
C NLRL    -- nodeElementOffsets
C LRL     -- nodeElementsArray
C IDPT    -- flag, a little confusing because combines node id and flag values
C XPT     -- x_out
C TPT     -- x_arrive_times
C MPT     -- element_track
C DT_INIT0-- target intial time step
C IDVE    -- form of velocity
C ID_DT   -- flag for initial time step choice (0 DT_INIT0, 1 TRY CFL=1)

      SUBROUTINE PT123
     I     (IDVE,IVERBOSE,MAXEQ,
     I      NNP,NEL,NNDE,NEQ,NPT,IBF,LU_O,
     I      ATOL,RTOL,SF,DN_SAFE,
     I      DT_INIT0,ID_DT,ID_RK,
     I      XG,IE,CL,NLRL,LRL,IB,
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
C ELEMENT DIAMETERS ARRAY
      DOUBLE PRECISION CL(NEL)
Cf2py double precision, intent (in) :: CL(NEL) 
C NODE - ELEMENTS IN NODE STAR LOOKUP (nodeElementOffsets,nodeElementsArray)
      INTEGER LRL(*),NLRL(*)
Cf2py integer, intent (in)  :: LRL(*),NLRL(*)
C FLAG ARRAY TO MARK NODES THAT ARE ON EXTERIOR BOUNDARY
      INTEGER IB(NNP)
Cf2py integer, intent (in)  :: IB(NNP)
C --TRACKING PARAMETERS--
      DOUBLE PRECISION ATOL,RTOL,SF,DN_SAFE,DT_INIT0
Cf2py double precision, intent (in) :: ATOL,RTOL,SF,DN_SAFE,DT_INIT0
C HOW TO PICK INITIAL TIME STEP ON EACH ELEMENT 0 --> CONSTANT, 1 --> CFL=1
      INTEGER ID_DT
Cf2py integer, intent (in) :: ID_DT
C TYPE OF RK TO USE, 45, 24,  
      INTEGER ID_RK
Cf2py integer, intent (in) :: ID_RK

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
CMWF TODO
C
C 
C -- LOCAL VARIABLES --
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
      DOUBLE PRECISION DN_S(MAXND),DN(MAXND)
C
      DOUBLE PRECISION XOUTB(MAXEQ),XOUTA(MAXEQ),XERR(MAXEQ)
      DOUBLE PRECISION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ)
      DOUBLE PRECISION AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ)
C 
      DOUBLE PRECISION XS(MAXEQ),XTEMP(MAXEQ)
      DOUBLE PRECISION DEQ,TS,DTS,T,TT,SDT,DT0
      DOUBLE PRECISION DIR,DL,DT_INIT
C     
      INTEGER I,IPT,NODE,NP,M,K,ID_RK0,IPROJ,ID_ETA
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
CMWF PROTEUS STILL ASSUMES NODE = NSPACE+1 FOR NOW
      NODE = NNDE
      DIR = 1.D0
      IF (IBF.EQ.-1) THEN
         DIR = -1.D0
      ENDIF
CMWF ALWAYS TIME INTERPOLATE FOR NOW
      ID_ETA = 1
      ID_RK0 = ID_RK
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
        ID_RK=ID_RK0
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
        
        IPROJ=0
CMWF DIFFERENT THAN PIERCE'S CONVENTION, IF NP==0 PROTEUS POINT IS IN THE INTERIOR, OTHERWISE IT'S A NODE
        IF(NP.EQ.0)GOTO 200

C
C ===== FOR THE CASE THAT MPT(IPT)=0, I.E., THE VERY FIRST
C     TRACKING: LOOP OVER ALL CONNECTED ELEMENTS
C
        IF (IVERBOSE.GE.4) THEN
C MWF DEBUG
           WRITE(LU_O,*)' ENTERING NODE TRACKING STEP IPT= ', IPT
        ENDIF
CMWF IB CONVENTION DIFFERENT FROM PEARCE, 1 --> ON BOUNDARY HERE, -1 ON BOUNDARY FOR PEARCE
        IF(IB(NP).EQ.1)IPROJ=1
  149   CONTINUE
        DO 150 I=NLRL(NP)+1,NLRL(NP+1)
          M=LRL(I)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,IPROJ,
     &         XG,IE,IB,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
CMWF NEW 
          DL=CL(M)
          IF(ID_DT.EQ.0)THEN
            DT_INIT=DT_INIT0
          ELSE
            CALL DT_ETRACK
     I          (MAXND,MAXEQ,NODE,NEQ, 
     I           DL,VT1W,VT2W,
     O           DT_INIT)
          ENDIF
          
C MWF REMOVE NPATH DIMENSIONS, CALL WITH MAXPT=NPT
          CALL ELTRAK123
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,M, ID_ETA,
     I         IPT, ID_RK,ID_RK0,T1,T2, DL, ATOL,RTOL,SF, DN_SAFE,
     I         XW,VT1W,VT2W,
     M         T,DT0,SDT,DT_INIT,XS,
     M         AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M         XOUTB,XOUTA,XERR,DN_S,DN,
     O         IDSDT,XPT,TPT,I1,I2,I3)
C          NPATH(IPT)=KPATH
C
          IF(IDSDT.EQ.-1)THEN
            MPT(IPT)=M
C
C 0 -- INTERIOR
            IDPT(IPT)=0
CMWF SKIP PATH OUTPUT
C            DO IEQ=1,NEQ
C              XPTW(IEQ,IPT)=XOUTB(IEQ)
C            ENDDO
            GOTO 600
          ENDIF
          IF(IDSDT.EQ.1)GOTO 250
  150   CONTINUE
C
C === SET ID_RK TO 1 TO CONTINUE PT
C
        IF(ID_RK.NE.1)THEN
          ID_RK=1
          GOTO 149
        ENDIF
C
C === NO TRACKING WAS PERFORMED, PRINT PARTICLE INFORMATION
C
C -1 -- EXITED DOMAIN
C -2 -- FAILED
C NOTE PEARCE CHECKS IF IB(NP) == 1 MEANING OPEN BOUNDARY, HERE WE SAY SUCCESS IF 
C JUST ON BOUNDARY AND DON'T DISTINGUISH OPEN OR CLOSED
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
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,IPROJ,
     &         XG,IE,IB,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)

C
C CONDUCT TRACKING WITHIN ELEMENT M
C
        DL=CL(M)
        IF(ID_DT.EQ.0)THEN
          DT_INIT=DT_INIT0
        ELSE
          CALL DT_ETRACK
     I        (MAXND,MAXEQ,NODE,NEQ, 
     I         DL,VT1W,VT2W,
     O         DT_INIT)
        ENDIF
C
        CALL ELTRAK123
     I      (MAXEQ,MAXND,NPT,NEQ,NODE,M,ID_ETA,
     I       IPT, ID_RK, ID_RK0, T1,T2, DL, ATOL,RTOL,SF, DN_SAFE,
     I       XW,VT1W,VT2W,
     M       T,DT0,SDT,DT_INIT,XS,
     M       AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M       XOUTB,XOUTA,XERR,DN_S,DN,
     O       IDSDT,XPT,TPT,I1,I2,I3)
        IF (IVERBOSE.GE.5) THEN
C MWF DEBUG
           WRITE(LU_O,*)' AFTER ELTRACK ELEMENT LOOP T= ',T,' SDT= ',
     &          SDT,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT,
     &          ' M= ',M
           WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
           DO K=1,NEQ
              WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &             XPT(K + (IPT-1)*MAXEQ)
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
        ID_RK=ID_RK0
        TT=T+SDT
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
          IPROJ=0
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1)IPROJ=1
  299     CONTINUE
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
     &                  NNP,NEL,NODE,NEQ,M,IDVE,DIR,IPROJ,
     &                  XG,IE,IB,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
                  DL=CL(M)
                  IF(ID_DT.EQ.0)THEN
                    DT_INIT=DT_INIT0
                  ELSE
                    CALL DT_ETRACK
     I                  (MAXND,MAXEQ,NODE,NEQ, 
     I                   DL,VT1W,VT2W,
     O                   DT_INIT)
                  ENDIF
C
                  CALL ELTRAK123
     I                (MAXEQ,MAXND,NPT,NEQ,NODE,M, ID_ETA,
     I                 IPT, ID_RK,ID_RK0, T1,T2, DL, ATOL,RTOL,SF,
     I                 DN_SAFE, XW,VT1W,VT2W,
     M                 T,DT0,SDT,DT_INIT,XS,
     M                 AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M                 XOUTB,XOUTA,XERR,DN_S,DN,
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
     &                       XPT(K + (IPT-1)*MAXEQ)
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
C === SET ID_RK TO 1 TO CONTINUE PT
C
          IF(ID_RK.NE.1)THEN
            ID_RK=1
            GOTO 299
          ENDIF
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
          IPROJ=0
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1)IPROJ=1
  399     CONTINUE
          DO J1=NLRL(N1)+1,NLRL(N1+1)          
            M=LRL(J1)
            DO J2=NLRL(N2)+1,NLRL(N2+1)          
              M2=LRL(J2)
              IF(M.EQ.M2)THEN
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
                 CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &                NNP,NEL,NODE,NEQ,M,IDVE,DIR,IPROJ,
     &                XG,IE,IB,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
                DL=CL(M)
                IF(ID_DT.EQ.0)THEN
                  DT_INIT=DT_INIT0
                ELSE
                  CALL DT_ETRACK
     I                (MAXND,MAXEQ,NODE,NEQ, 
     I                 DL,VT1W,VT2W,
     O                 DT_INIT)
                ENDIF
C
                CALL ELTRAK123
     I              (MAXEQ,MAXND,NPT,NEQ,NODE,M, ID_ETA,
     I               IPT, ID_RK, ID_RK0, T1,T2, DL, ATOL,RTOL,SF,
     I               DN_SAFE, XW,VT1W,VT2W,
     M               T,DT0,SDT,DT_INIT,XS,
     M               AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M               XOUTB,XOUTA,XERR,DN_S,DN,
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
     &                     XPT(K + (IPT-1)*MAXEQ)
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
C === SET ID_RK TO 1 TO CONTINUE PT
C
          IF(ID_RK.NE.1)THEN
            ID_RK=1
            GOTO 399
          ENDIF
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
        IPROJ=0
        IF(IB(N1).EQ.1)IPROJ=1
  499   CONTINUE
        DO J1=NLRL(N1)+1,NLRL(N1+1)          
          M=LRL(J1)
C ELENOD IN GENERAL SELECTS THE NUMBER OF NODES PER ELEMENT USING IE
C FOR NOW WE KNOW THIS AS INPUT
C
          CALL EL_VEL_PREP(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,IPROJ,
     &         XG,IE,IB,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)
C
C CONDUCT TRACKING WITHIN ELEMENT M
C
          DL=CL(M)
          IF(ID_DT.EQ.0)THEN
            DT_INIT=DT_INIT0
          ELSE
            CALL DT_ETRACK
     I          (MAXND,MAXEQ,NODE,NEQ, 
     I           DL,VT1W,VT2W,
     O           DT_INIT)
          ENDIF
C
          CALL ELTRAK123
     I        (MAXEQ,MAXND,NPT,NEQ,NODE,M, ID_ETA,
     I         IPT, ID_RK, ID_RK0, T1,T2, DL, ATOL,RTOL,SF, DN_SAFE,
     I         XW,VT1W,VT2W,
     M         T,DT0,SDT,DT_INIT,XS,
     M         AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M         XOUTB,XOUTA,XERR,DN_S,DN,
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
                WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &               XPT(K + (IPT-1)*MAXEQ)
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
C === SET ID_RK TO 1 TO CONTINUE PT
C
        IF(ID_RK.NE.1)THEN
          ID_RK=1
          GOTO 499
        ENDIF
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
     I    (MAXEQ,MAXND,MAXPT,NEQ,NODE,M, ID_ETA,
     I     IPT, ID_RK, ID_RK0, T1,T2,DL, ATOL,RTOL,SF, DN_SAFE,
     I     XW,VT1W,VT2W,
     M     T,DT0,SDT,DT_INIT, XS, 
     M     AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     M     XOUTB,XOUTA,XERR,DN_S,DN,
     O     IDSDT,XPT,TPT,I1,I2,I3)
C 
C 07/14/2010 (HPC)
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
      DIMENSION XOUTB(MAXEQ),XOUTA(MAXEQ),XERR(MAXEQ)
      DIMENSION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ),
     >          AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ),
     >          XTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND)
      DIMENSION ICHECK(8)
      DIMENSION XI(3),DI(3),XI_S(3),DI_S(3)
      DIMENSION IXI(3),IDI(4),IXI_S(3),IDI_S(4)
CMWF LOCAL VARIABLES
      INTEGER IDEBUG,II,KK,ISTEPOUTSIDE
C
C =================== START PT USING ADAPTIVE RK ====================
C
C TRACKING IN ELEMENT M
C
      IDEBUG=0
      ICOUNT=0
      IDSDT=1
      DT=DT_INIT*DSIGN(1.D0,DT0)
      IF (IDEBUG.GT.1) THEN
CMWF DEBUG
         WRITE(6,*)'ENTERING ELTRAK T= ',T,'DT0= ',DT0,' DT= ',DT,
     &     ' SDT= ',SDT,' T1= ',T1,' T2= ',T2,' ID_RK= ',ID_RK,
     &     ' ID_RK0= ',ID_RK0,' DN_SAFE= ',DN_SAFE
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
      ENDIF
CMWF DEBUG

      IF(ID_RK.NE.ID_RK0)THEN
        DT1=1.0E-3*DT
        DT2=1.0E1*DN_SAFE*DL*DSIGN(1.D0,DT0)
CMWF        DT=DMAX1(DT1,DT2)
        DT = DT1
        IF (DABS(DT2).GT.DABS(DT1)) THEN 
           DT = DT2
        ENDIF
      ENDIF
C
  100 CONTINUE
      ICOUNT=ICOUNT+1
      IF(ID_RK.EQ.45)THEN
        CALL RKCK_PT
     I      (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, T,DT,T1,T2,
     I       XW,VT1W,VT2W,XS, ID_ETA,
     M       AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     O       IDSDT,IREDUCE,I1,I2,I3,
     O       XOUTB,XOUTA,XERR,DN_S,XI_S,DI_S,IXI_S,IDI_S)
        P1=0.2E0
        P2=0.25E0
      ELSEIF(ID_RK.LE.24)THEN
        CALL RK24_PT
     I      (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, T,DT, T1,T2,
     I       XW,VT1W,VT2W, XS, ID_ETA, ID_RK,
     M       AK1,AK2,AK3,AK4,XTEMP,
     O       IDSDT,IREDUCE,I1,I2,I3,
     O       XOUTB,XOUTA,XERR,DN_S,XI_S,DI_S,IXI_S,IDI_S)
        P1=0.35E0
        P2=0.40E0
      ENDIF
C
      IF (IDEBUG.GT.1) THEN
CMWF DEBUG
         WRITE(6,*)'ELTRAK AFTER RKCK_PT T= ',T,' DT= ',DT,
     &     ' IDSDT= ',IDSDT,' IREDUCE= ',IREDUCE
         WRITE(6,*)'XOUTB= ',(XOUTB(I),I=1,NEQ)
         WRITE(6,*)'DN_S= ',(DN_S(I),I=1,NODE)
      ENDIF
CMWF DEBUG
      IF(IDSDT.EQ.0)RETURN
      IF(IREDUCE.EQ.1)GOTO 100
C
C CHECK ERROR IN ALL THREE DIRECTIONS 
C
      RATIO=0.0D0
      DO I=1,NEQ
        XERRABS=DABS(XERR(I))
        XOUTBABS=DABS(XOUTB(I)-XS(I))
        XOUTAABS=DABS(XOUTA(I)-XS(I))
        XOUTMAX=DMAX1(XOUTBABS,XOUTAABS)
        RATIOI=XERRABS/(RTOL*XOUTMAX+ATOL)
        RATIO=DMAX1(RATIO,RATIOI)
      ENDDO
C
C CASE 1: WHEN RATIO IS GREATER THAN 1
C      ==> DESIRED ACCURACY IS NOT REACHED
C      ==> TIMESTEP NEEDS TO BE REDUCED
C
      IF(RATIO.GT.1.0E0)THEN
        RR=1.0E0/RATIO
        DTT=DT*SF*((RR)**P2) 
        IF (IDEBUG.GT.2) THEN
CMWF DEBUG
            WRITE(6,*)'ELTRAK AFTER ERROR FAILURE T= ',T,' RATIO= ',
     &           RATIO,' DT= ',DT,' DTT= ',DT*SF*((RR)**(0.25E0)) 
        ENDIF
CMWF DEBUG
        DT=DTT
        GOTO 100
C
C CASE 2: WHEN RATIO IS LESS THAN OR EQUAL TO 1
C      ==> DESIRED ACCURACY IS REACHED
C      ==> TIMESTEP CAN BE INCREASED
C        
      ELSEIF(RATIO.LE.1.0E0)THEN
CMWF DEBUG
         IF(IDEBUG.GT.2) THEN
            WRITE(6,*)'ELTRAK AFTER ERROR SUCCESS T= ',T,' RATIO= ',
     &           RATIO,' DT= ',DT,' DT0= ',DT0
         ENDIF
CMWF DEBUG
C DT CANNOT BE GREATER THAN DT0
C
CMWF ORIG        IF(DT.GT.DT0)THEN
        IF(DABS(DT).GT.DABS(DT0)) THEN
          DT=DT0
          GOTO 100
        ENDIF
C
C === EXAMINE THE COMPUTED ENDING LOCATION
C
        PHI=0.0D0
        CALL INTRP123
     I      (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XOUTB, XW,
     O       DN,IADJUST,XI,DI,IXI,IDI)
C
C < NOTE > ADJUST XOUTB WHEN NECESSARY (I.E., IADJUST = 1)
C
        IF (IDEBUG.GT.2) THEN 
           WRITE(6,*)'AFTER SUCCESS INTRP123 IADJUST= ',IADJUST
           WRITE(6,*)'DN= ',(DN(I),I=1,NODE)
        ENDIF
        IF(IADJUST.EQ.1)THEN
          DO I=1,NEQ
            XOUTB(I)=0.0E0
            DO J=1,NODE
              XOUTB(I)=XOUTB(I)+DN(J)*XW(I,J)
            ENDDO
          ENDDO
        ENDIF                 
C
C === COMPUTE PHI
C
        CALL PHI_COMP
     I      (MAXND,NODE,NEQ, 
     I       DN_SAFE, 
     I       DN_S,DN,XI,XI_S,DI,DI_S, IXI_S,IDI_S,
     O       IDSDT,I1,I2,I3,PHI)
CMWF
        IF(IDEBUG.GT.2) THEN
           WRITE(6,*)'AFTER ERROR SUCCESS PHI_COMP IDSDT= ',IDSDT,
     &          ' PHI= ',PHI
        ENDIF

        IF(IDSDT.EQ.0)RETURN
C
C === IF PHI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C
        IF(PHI.GT.1.0E0)THEN
          DTT=DT/PHI
          DT=DTT
          GOTO 100
        ENDIF
C
C B. WHEN THE ENDING LOCATION IS EITHER WITHIN THE ELEMENT OR 
C    ON THE ELEMENT BOUDNARY
C
        T=T+DT
        SDT=SDT-DT
        DT0=DT0-DT
        IF(RATIO.LT.1.0E-6)THEN
          DTT=DT
        ELSE
          RR=1.0E0/RATIO
          DTT=DT*SF*((RR)**P1)
        ENDIF
CMWF DEBUG
        IF(IDEBUG.GT.2) THEN
           WRITE(6,*)'AFTER FULL STEP T= ',T,' SDT= ',SDT,
     &          ' DT0= ',DT0
           WRITE(6,*)'XOUTB= ',(XOUTB(I),I=1,NEQ)
        ENDIF
C
C ... IF THE LOCATION CHANGE OF PARTICLE IS NEGLIGIBLE
C     ===> SET IDSDT = 0
C
c        DIFF=0.0E0
c        DO I=1,NEQ
c          DD=XOUTB(I)-XS(I)
c          DIFF=DIFF+DD*DD
c        ENDDO
c        IF(DSQRT(DIFF).LT.ATOL)IDSDT=0
C
C ... STORE TRACKING LOCATIONS AT THE SPECIFIED FREQUENCY
C
C       print *,'t, id_rk =',t,id_rk
C
CMWF remove storage
C        NTPATH(IPT)=NTPATH(IPT)+1
C        CALL PT_STORE
C     I    (MAXEQ,MAXPATH,MAXPT, NEQ, ID_RK,
C     I     IPT,T,DT,DT_OUTPUT,T_OUTPUT,
C     I     XOUTB,XS,
C     M     KPATH,XPT,TPT,ID_RKPT)
C
C =======================================
C FOR DEBUGGING ONLY
C
c        if(kpath.ge.515)then
c           print *,'kpath=',kpath
c        endif
C
C =======================================
C
C ... UPDATE INFORMATION FOR THE SUCCESSIVE PT
C
        DO I=1,NEQ
          XS(I)=XOUTB(I)
CMWF CHANGED XPT TO TAKE CARE OF INPUT OUTPUT OF POINTS
          XPT(I + MAXEQ*(IPT-1))=XS(I)
          
        ENDDO  
CMWF ADDED
        TPT(IPT) = T
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
        CALL EB_CHECK
     I      (NEQ,NODE, IDI,IXI, 
     O       I1,I2,I3)
        IF(I1.NE.0)RETURN
C
C UPDATE THE TRACKING TIMESTEP (I.E., DT) FOR THE SUCCESSIVE TRACKING
C WITHIN THIS SAME ELEMENT
C
        IF(IDSDT.EQ.-1)RETURN
        IF(ID_RK.EQ.24 .OR. ID_RK.EQ.45)THEN
          DT=DT0
CMWF ORIG          DT=DMIN1(DT,DTT)
          IF(DABS(DTT).LT.DABS(DT)) THEN
             DT=DTT
          ENDIF
        ELSE
          DT=DT_INIT*DSIGN(1.D0,DT0)
        ENDIF
        ID_RK=ID_RK0
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
      SUBROUTINE RKCK_PT
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL,T,DT, T1,T2,
     I     XW,VT1W,VT2W, XS,ID_ETA,
     M     AK1,AK2,AK3,AK4,AK5,AK6,XTEMP,
     O     IDSDT,IREDUCE,I1,I2,I3,
     O     XOUTB,XOUTA,XERR,DN_S,XI_S,DI_S,IXI_S,IDI_S)
C 
C 07/14/2010 (HPC)
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
C   XOUTB  = ESTIMATE FROM 5TH-ORDER RK
C   XOUTA  = ESTIMATE FROM 4TH-ORDER RK
C   XERR   = ERROR ESTIMATE BETWEEN 4TH- AND 5TH-ORDER RK
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XS(MAXEQ)
      DIMENSION XOUTB(MAXEQ),XOUTA(MAXEQ),XERR(MAXEQ)
      DIMENSION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ),
     >          AK4(MAXEQ),AK5(MAXEQ),AK6(MAXEQ),
     >          XTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND)
      DIMENSION XI(3),DI(3),XI_S(3),DI_S(3)
      DIMENSION IXI(3),IDI(4),IXI_S(3),IDI_S(4)
CMWF FOR DEBUGGING
      INTEGER IDEBUG
C
      IDEBUG = 0
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
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK1,DN_S,XI_S,DI_S,IXI_S,IDI_S)
C
C === EXAMINE THE COMPUTED ENDING LOCATION TO DETERMINE WHETHER THIS IS THE 
C     ELEMENT THE PARTICLE WILL TRAVEL WITHIN IT
C
      DO I=1,NEQ
        XOUTA(I)=XS(I)+DT*AK1(I)
      ENDDO
      IREDUCE=0
      PHI=0.0D0
      CALL INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XOUTA, XW,
     O     DN,IADJUST,XI,DI,IXI,IDI)
      CALL PHI_COMP
     I    (MAXND,NODE,NEQ, 
     I     DN_SAFE,  
     I     DN_S,DN,XI,XI_S,DI,DI_S, IXI_S,IDI_S,
     O     IDSDT,I1,I2,I3,PHI)
C
CMWF DEBUG
      IF(IDEBUG.GT.1) THEN
         WRITE(6,*)'RKCK_PT T= ',T,' TT= ',TT,' DT= ',DT,
     &        ' IDSDT= ',IDSDT,' I1= ',I1,' I2= ',I2,' I3= ',I3,
     &        ' PHI= ',PHI
         WRITE(6,*)'AK1= ',(AK1(I),I=1,NEQ)
         WRITE(6,*)'XOUTA= ',(XOUTA(I),I=1,NEQ)
         WRITE(6,*)'DN= ',(DN(I),I=1,NODE)
      ENDIF
      IF(IDSDT.EQ.0)RETURN
C
C === IF PHI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C
      IF(PHI.GT.1.0E0)THEN
        DTT=DT/PHI
        DT=DTT
        IREDUCE=1
        IF(IDEBUG.GT.1) THEN
         WRITE(6,*)'RKCK_PT REDUCING, PHI= ',PHI,
     &        ' DTT= ',DTT,' IREDUCE= ',IREDUCE

      ENDIF

        RETURN
      ENDIF
C
C STEP 1:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*B21*AK1(I)
      ENDDO
      TT=T+A2*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK2,DN,XI,DI,IXI,IDI)
C
C STEP 2:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B31*AK1(I)+B32*AK2(I))
      ENDDO
      TT=T+A3*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK3,DN,XI,DI,IXI,IDI)
C
C STEP 3:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B41*AK1(I)+B42*AK2(I)+B43*AK3(I))
      ENDDO
      TT=T+A4*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK4,DN,XI,DI,IXI,IDI)
C
C STEP 4:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B51*AK1(I)+B52*AK2(I)+B53*AK3(I)+
     >                     B54*AK4(I))
      ENDDO
      TT=T+A5*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK5,DN,XI,DI,IXI,IDI)
C
C STEP 5:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*(B61*AK1(I)+B62*AK2(I)+B63*AK3(I)+
     >                     B64*AK4(I)+B65*AK5(I))
      ENDDO
      TT=T+A6*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK6,DN,XI,DI,IXI,IDI)
C
C ESTIMATE ERROR USING 5TH- AND 4TH-ORDER RK
C
      DO I=1,NEQ
        XOUTB(I)=XS(I)+DT*(C1*AK1(I)+C2*AK2(I)+C3*AK3(I)+
     >                     C4*AK4(I)+C5*AK5(I)+C6*AK6(I))
        XOUTA(I)=XS(I)+DT*(D1*AK1(I)+D2*AK2(I)+D3*AK3(I)+
     >                     D4*AK4(I)+D5*AK5(I)+D6*AK6(I))
        XERR(I)=DT*(E1*AK1(I)+E2*AK2(I)+E3*AK3(I)+
     >              E4*AK4(I)+E5*AK5(I)+E6*AK6(I))
      ENDDO
C
C ===== RETURN TO THE CALLING ROUTINE
C
  999 CONTINUE
      RETURN
      END
C
C 
C
      SUBROUTINE RK24_PT
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, T,DT, T1,T2,
     I     XW,VT1W,VT2W, XS, ID_ETA, ID_RK,
     M     AK1,AK2,AK3,AK4,XTEMP,
     O     IDSDT,IREDUCE,I1,I2,I3,
     O     XOUTB,XOUTA,XERR,DN_S,XI_S,DI_S,IXI_S,IDI_S)
C
C 07/14/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C   ESTIMATE ERRORS BETWEEN 2ND- AND 4TH-ORDER RUNGE-KUTTA (RK) METHOD 
C   CONCERNING PARTICLE TRACKING (PK)
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
C   XOUTB  = ESTIMATE FROM 2ND-ORDER RK
C   XOUTA  = ESTIMATE FROM 4TH-ORDER RK
C   XERR   = ERROR ESTIMATE BETWEEN 4TH- AND 5TH-ORDER RK
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XS(MAXEQ)
      DIMENSION XOUTB(MAXEQ),XOUTA(MAXEQ),XERR(MAXEQ)
      DIMENSION AK1(MAXEQ),AK2(MAXEQ),AK3(MAXEQ),
     >          AK4(MAXEQ),XTEMP(MAXEQ)
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
      DIMENSION DN_S(MAXND),DN(MAXND)
      DIMENSION XI(3),DI(3),XI_S(3),DI_S(3)
      DIMENSION IXI(3),IDI(4),IXI_S(3),IDI_S(4)
C
C ===== DEFINE THE PARAMETERS USED FOR 2ND- AND 4-TH
C       ORDER RK
C
      D1  = 1.0E0/6.0E0
      D2  = 1.0E0/3.0E0
      D3  = D2
      D4  = D1
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
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK1,DN_S,XI_S,DI_S,IXI_S,IDI_S)
C
C === EXAMINE THE COMPUTED ENDING LOCATION TO DETERMINE WHETHER THIS IS THE 
C     ELEMENT THE PARTICLE WILL TRAVEL WITHIN IT
C
      DO I=1,NEQ
        XOUTA(I)=XS(I)+DT*AK1(I)
      ENDDO
      IREDUCE=0
      PHI=0.0D0
      CALL INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XOUTA, XW,
     O     DN,IADJUST,XI,DI,IXI,IDI)
      CALL PHI_COMP
     I    (MAXND,NODE,NEQ, 
     I     DN_SAFE,  
     I     DN_S,DN,XI,XI_S,DI,DI_S, IXI_S,IDI_S,
     O     IDSDT,I1,I2,I3,PHI)
C
      IF(IDSDT.EQ.0)RETURN
C
C === IF PHI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C
      IF(PHI.GT.1.0E0)THEN
        DTT=DT/PHI
        DT=DTT
        IREDUCE=1
        RETURN
      ENDIF
C
C WHEN 1ST-ORDER RK IS DESIRED
C
      IF(ID_RK.EQ.1)THEN
        DO I=1,NEQ
          XOUTB(I)=XOUTA(I)
          XERR(I)=0.0E0
        ENDDO
        RETURN
      ENDIF
C
C STEP 1:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*0.5E0*AK1(I)
      ENDDO
      TT=T+0.5E0*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK2,DN,XI,DI,IXI,IDI)
C
C WHEN 2ND-ORDER RK IS DESIRED
C
      IF(ID_RK.EQ.2)THEN
        DO I=1,NEQ
          XOUTA(I)=XS(I)+DT*AK2(I)
          XOUTB(I)=XOUTA(I)
          XERR(I)=0.0E0
        ENDDO
        RETURN
      ENDIF
C
C STEP 2:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*0.5E0*AK2(I)
      ENDDO
      TT=T+0.5E0*DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK3,DN,XI,DI,IXI,IDI)
C
C STEP 3:
C
      DO I=1,NEQ
        XTEMP(I)=XS(I)+DT*AK3(I)
      ENDDO
      TT=T+DT
      CALL VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK4,DN,XI,DI,IXI,IDI)
C
C WHEN 4ND-ORDER RK IS DESIRED
C
      IF(ID_RK.EQ.4)THEN
        DO I=1,NEQ
          XOUTA(I)=XS(I)+DT*(D1*AK1(I)+D2*AK2(I)+
     >                     D3*AK3(I)+D4*AK4(I))
          XOUTB(I)=XOUTA(I)
          XERR(I)=0.0E0
        ENDDO
        RETURN
      ENDIF
C
C ESTIMATE ERROR USING 2ND- AND 4TH-ORDER RK
C
      DO I=1,NEQ
        XOUTA(I)=XS(I)+DT*AK2(I)
        XOUTB(I)=XS(I)+DT*(D1*AK1(I)+D2*AK2(I)+
     >                     D3*AK3(I)+D4*AK4(I))
        XERR(I)=XOUTA(I)-XOUTB(I)
      ENDDO
C
C ===== RETURN TO THE CALLING ROUTINE
C
  999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE VEL123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP,TT, T1,T2,
     I     XW,VT1W,VT2W, ID_ETA,
     O     AK,DN,XI,DI,IXI,IDI)
C 
C 06/22/2010 (HPC)
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
      DIMENSION XI(3),DI(3),IXI(3),IDI(4)
C
C
C === COMPUTE THE INTERPOLATION FUNCTIONAL VALUES
C
C IN SPACE:
C
      CALL INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP, XW,
     O     DN,IADJUST,XI,DI,IXI,IDI)
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
C      WRITE(6,*)' VEL123 T1= ',T1, 'T2= ',T2,' TT= ',TT,' ETA= ',ETA
C      DO I=1,NEQ
C         WRITE(6,*)' VT1(',I,')= ',VT1(I),' VT2(',I,')= ',VT2(I),' AK(',
C     &        I,')= ',AK(I)
C      ENDDO
C
C  999 CONTINUE
      RETURN
      END 
C
C
C
      SUBROUTINE INTRP123
     I    (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ATOL,DL, XTEMP, XW,
     O     DN,IADJUST,XI,DI,IXI,IDI)
C
C 06/22/2010 (HPC) 
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
C
C < NOTE >
C NEQ = 1    ==> 1-D LINE ELEMENT
C NEQ = 2
C   NODE = 3 ==> 2-D TRIANGULAR ELEMENT
C   NODE = 4 ==> 2-D QUADRILATERAL ELEMENT
C NEQ = 3  
C   NODE = 4 ==> 3-D TETRAHEDRAL ELEMENT
C   NODE = 6 ==> 3-D TRIANGULAR PRISM ELEMENT
C   NODE = 8 ==> 3-D HEXAHEDRAL ELEMENT
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
      DIMENSION DN(MAXND),XI(3),DI(3),IXI(3),IDI(4)
      DIMENSION A(4),B(4),C(4),D(4)
      DIMENSION K21(3,3),K31(4,4),K32(4,6)
      DATA K21 /1,2,3, 2,3,1, 3,1,2/
      DATA K31 /1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3/
      DATA K32 /1,2,3,4, 1,3,2,4, 1,4,2,3, 2,3,1,4, 2,4,1,3, 3,4,1,2/    
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
        XSI=(XTEMP(1)-XW(1,1))/(XW(1,2)-XW(1,1))
        DN(1)=1.0E0-XSI
        DN(2)=XSI
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)
C
C
C ===== FOR THE CASE OF A 2-D TRIANGULAR ELEMENT
C
      ELSEIF(NEQ.EQ.2 .AND. NODE.EQ.3)THEN
        X12=XW(1,1)-XW(1,2)
        X23=XW(1,2)-XW(1,3)
        X31=XW(1,3)-XW(1,1)  
        Y12=XW(2,1)-XW(2,2)
        Y23=XW(2,2)-XW(2,3)
        Y31=XW(2,3)-XW(2,1)
        DJAC=XW(1,1)*Y23+XW(1,2)*Y31+XW(1,3)*Y12
C
        DN(1)=(Y23*XTEMP(1)-X23*XTEMP(2)+XW(1,2)*XW(2,3)-
     >                                   XW(1,3)*XW(2,2))/DJAC
        DN(2)=(Y31*XTEMP(1)-X31*XTEMP(2)+XW(1,3)*XW(2,1)-
     >                                   XW(1,1)*XW(2,3))/DJAC
        DN(3)=1.0E0-DN(1)-DN(2)
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)         
C
C
C ===== FOR THE CASE OF A 2-D QUADRILATERL ELEMENT
C
      ELSEIF(NEQ.EQ.2 .AND. NODE.EQ.4)THEN
        CALL XSI_2 
     I      (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O       XSI,ETA)
        XI(1)=XSI
        XI(2)=ETA
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)
C
C
C ===== FOR THE CASE OF A 3-D TETRAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.4)THEN
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
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)
C
C
C ===== FOR THE CASE OF A 3-D TRIANGULAR PRISM ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.6)THEN
        CALL XSI_3P
     I      (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O       XSI,DL1,DL2,DL3)
        XI(1)=XSI
        DI(1)=DL1
        DI(2)=DL2
        DI(3)=DL3
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)
C
C
C ===== FOR THE CASE OF A 3-D HEXAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.8)THEN
        CALL XSI_3
     I      (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O       XSI,ETA,ZTA)
        XI(1)=XSI
        XI(2)=ETA
        XI(3)=ZTA
C
        CALL ADJUST123
     I      (MAXND,NEQ,NODE,DN_SAFE,
     M       DN,XI,DI,
     O       IADJUST,IXI,IDI)
C
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE ADJUST123
     I    (MAXND,NEQ,NODE,DN_SAFE,
     M     DN,XI,DI,
     O     IADJUST,IXI,IDI)
C
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C (1) DETERMINE IADJUST, IXI, AND IDI FOR NECESSARY LOCATION SHIFT WHEN
C     THE PARTICLE IS SUFFICIENTLY CLOSE TO THE ELEMENT BOUNDARY
C (2) ADJUST DN, XI, AND DI AS NECESSARY 
C
C < NOTE >
C NEQ = 1    ==> 1-D LINE ELEMENT
C NEQ = 2
C   NODE = 3 ==> 2-D TRIANGULAR ELEMENT
C   NODE = 4 ==> 2-D QUADRILATERAL ELEMENT
C NEQ = 3  
C   NODE = 4 ==> 3-D TETRAHEDRAL ELEMENT
C   NODE = 6 ==> 3-D TRIANGULAR PRISM ELEMENT
C   NODE = 8 ==> 3-D HEXAHEDRAL ELEMENT
C
C   DI(I) = NATURAL COORDINATE ASSOCAITED THE I-TH ELEMENT NODE [0,1]
C   IDI(I) = -1, WHEN XI(I) = 0
C             1, WHEN XI(I) = 1
C             0, OTHERWISE
C   XI(I) = THE I-TH LOCAL COORDINATE [-1,1]
C   IXI(I) = -1, WHEN XI(I) = -1
C             1, WHEN XI(I) = 1
C             0, OTHERWISE
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DN(MAXND),XI(3),DI(3),IXI(3),IDI(4)
      DIMENSION K21(3,3),K31(4,4),K32(4,6)
      DATA K21 /1,2,3, 2,3,1, 3,1,2/
      DATA K31 /1,2,3,4, 2,3,4,1, 3,4,1,2, 4,1,2,3/
      DATA K32 /1,2,3,4, 1,3,2,4, 1,4,2,3, 2,3,1,4, 2,4,1,3, 3,4,1,2/    
C
      IADJUST=0
      DN_SAFE2=2.0E0*DN_SAFE
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
        IF(DABS(DN(1)).LE.DN_SAFE)THEN
          IADJUST=1
          DN(1)=0.0E0
          DN(2)=1.0E0
          IDI(1)=-1
          IDI(2)=1
        ELSEIF(DABS(DN(1)-1.0E0).LE.DN_SAFE)THEN
          IADJUST=1
          DN(1)=1.0E0
          DN(2)=0.0E0
          IDI(1)=1
          IDI(2)=-1
        ELSE
          IDI(1)=0
          IDI(2)=0
        ENDIF
C
C
C ===== FOR THE CASE OF A 2-D TRIANGULAR ELEMENT
C
      ELSEIF(NEQ.EQ.2 .AND. NODE.EQ.3)THEN
        ICOUNT=0
        DO I=1,3
          IDI(I)=0
          IF(DABS(DN(I)).LE.DN_SAFE)THEN
            DN(I)=0.0E0
            IDI(I)=-1
            ICOUNT=ICOUNT+1
          ELSEIF(DABS(DN(I)-1.0E0).LE.DN_SAFE)THEN
            DN(I)=1.0E0
            IDI(I)=1
            ICOUNT=ICOUNT+1
          ENDIF
        ENDDO
        IF(ICOUNT.EQ.3)THEN
          IADJUST=1
        ELSEIF(ICOUNT.EQ.2)THEN
          IADJUST=1
          DO I=1,3
            K1=K21(1,I)
            K2=K21(2,I)
            K3=K21(3,I)
            IF(IDI(K1).NE.-1 .AND. IDI(K1).NE.1)THEN
              DN(K1)=1.0E0-DN(K2)-DN(K3)
              IDI(K1)=1
              IF(DABS(DN(K1)).LE.DN_SAFE)IDI(K1)=-1
              RETURN
            ENDIF
          ENDDO
        ELSEIF(ICOUNT.EQ.1)THEN
          IADJUST=1
          DO I=1,3
            K1=K21(1,I)
            K2=K21(2,I)
            K3=K21(3,I)
            IF(IDI(K1).EQ.-1 .OR. IDI(K1).EQ.1)THEN
              DIFF=1.0E0-DN(K1)-DN(K2)-DN(K3)
              DN(K2)=DN(K2)+DIFF
              IF(DABS(DN(K2)).LE.DN_SAFE)THEN
                DN(K2)=0.0E0
                IDI(K2)=-1
                DN(K3)=1.0E0-DN(K1)-DN(K2)
                IDI(K3)=1
                IF(DABS(DN(K3)).LE.DN_SAFE)IDI(K3)=-1
                RETURN
              ELSEIF(DABS(DN(K2)-1.0E0).LE.DN_SAFE)THEN
                DN(K2)=1.0E0
                IDI(K2)=1
                DN(K3)=0.0E0
                IDI(K3)=-1
                RETURN
              ENDIF
            ENDIF
          ENDDO
        ENDIF          
C
C
C ===== FOR THE CASE OF A 2-D QUADRILATERL ELEMENT
C
      ELSEIF(NEQ.EQ.2 .AND. NODE.EQ.4)THEN
        IXI(1)=0
        IF(DABS(XI(1)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=1.0E0
          IXI(1)=1
        ELSEIF(DABS(XI(1)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=-1.0E0
          IXI(1)=-1
        ENDIF
        IXI(2)=0
        IF(DABS(XI(2)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(2)=1.0E0
          IXI(2)=1
        ELSEIF(DABS(XI(2)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(2)=-1.0E0
          IXI(2)=-1
        ENDIF
        DN(1)=0.25E0*(1.0E0-XI(1))*(1.0E0-XI(2))
        DN(2)=0.25E0*(1.0E0+XI(1))*(1.0E0-XI(2))
        DN(3)=0.25E0*(1.0E0+XI(1))*(1.0E0+XI(2))
        DN(4)=0.25E0*(1.0E0-XI(1))*(1.0E0+XI(2))
C
C
C ===== FOR THE CASE OF A 3-D TETRAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.4)THEN
        ICOUNT=0
        DO I=1,4
          IDI(I)=0
          IF(DABS(DN(I)).LE.DN_SAFE)THEN
            DN(I)=0.0E0
            IDI(I)=-1
            ICOUNT=ICOUNT+1
          ELSEIF(DABS(DN(I)-1.0E0).LE.DN_SAFE)THEN
            DN(I)=1.0E0
            IDI(I)=1
            ICOUNT=ICOUNT+1
          ENDIF
        ENDDO
C
        IF(ICOUNT.EQ.4)THEN
          IADJUST=1
        ELSEIF(ICOUNT.EQ.3)THEN
          IADJUST=1
          DO I=1,4
            K1=K31(1,I)
            K2=K31(2,I)
            K3=K31(3,I)
            K4=K31(4,I)
            IF(IDI(K1).NE.-1 .AND. IDI(K1).NE.1)THEN
              DN(K1)=1.0E0-DN(K2)-DN(K3)-DN(K4)
              IDI(K1)=1
              IF(DABS(DN(K1)).LE.DN_SAFE)IDI(K1)=-1
              RETURN
            ENDIF
          ENDDO
        ELSEIF(ICOUNT.EQ.2)THEN
          IADJUST=1
          DO I=1,6
            K1=K32(1,I)
            K2=K32(2,I)
            K3=K32(3,I)
            K4=K32(4,I)
            IF(IDI(K1).EQ.-1 .AND. IDI(K1).EQ.1)THEN
              IF(IDI(K2).EQ.-1 .AND. IDI(K2).EQ.1)THEN
                DIFF=1.0E0-DN(K1)-DN(K2)-DN(K3)-DN(K4)
                DN(K3)=DN(K3)+DIFF
                IF(DABS(DN(K3)).LE.DN_SAFE)THEN
                  DN(K3)=0.0E0
                  IDI(K3)=-1
                  DN(K4)=1.0E0-DN(K1)-DN(K2)-DN(K3)
                  IDI(K4)=1
                  IF(DABS(DN(K4)).LE.DN_SAFE)IDI(K4)=-1
                  RETURN
                ELSEIF(DABS(DN(K3)-1.0E0).LE.DN_SAFE)THEN
                  DN(K3)=1.0E0
                  IDI(K3)=1
                  DN(K4)=0.0E0
                  IDI(K4)=-1
                  RETURN
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ELSEIF(ICOUNT.EQ.1)THEN
          IADJUST=1
          DO I=1,4
            K1=K31(1,I)
            K2=K31(2,I)
            K3=K31(3,I)
            K4=K31(4,I)
            IF(IDI(K1).EQ.-1 .OR. IDI(K1).EQ.1)THEN
              DIFF=1.0E0-DN(K1)-DN(K2)-DN(K3)-DN(K4)
              DN(K2)=DN(K2)+DIFF
              ICHECK=0
              IF(DABS(DN(K2)).LE.DN_SAFE)THEN
                DN(K2)=0.0E0
                IDI(K2)=-1
                ICHECK=1
              ELSEIF(DABS(DN(K2)-1.0E0).LE.DN_SAFE)THEN
                DN(K2)=1.0E0
                IDI(K2)=1
                ICHECK=1
              ENDIF
              IF(ICHECK.EQ.1)THEN
                DIFF=1.0E0-DN(K1)-DN(K2)-DN(K3)-DN(K4)
                DN(K3)=DN(K3)+DIFF
                IF(DABS(DN(K3)).LE.DN_SAFE)THEN
                  DN(K3)=0.0E0
                  IDI(K3)=-1
                  DN(K4)=1.0E0-DN(K1)-DN(K2)-DN(K3)
                  IDI(K4)=1
                  IF(DABS(DN(K4)).LE.DN_SAFE)IDI(K4)=-1
                  RETURN
                ELSEIF(DABS(DN(K3)-1.0E0).LE.DN_SAFE)THEN
                  DN(K3)=1.0E0
                  IDI(K3)=1
                  DN(K4)=0.0E0
                  IDI(K4)=-1
                  RETURN
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C
C ===== FOR THE CASE OF A 3-D TRIANGULAR PRISM ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.6)THEN
        IXI(1)=0
        IF(DABS(XI(1)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=1.0E0
          IXI(1)=1
        ELSEIF(DABS(XI(1)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=-1.0E0
          IXI(1)=-1
        ENDIF
C
        ICOUNT=0
        DO I=1,3
          IDI(I)=0
          IF(DABS(DI(I)).LE.DN_SAFE)THEN
            DI(I)=0.0E0
            IDI(I)=-1
            ICOUNT=ICOUNT+1
          ELSEIF(DABS(DI(I)-1.0E0).LE.DN_SAFE)THEN
            DI(I)=1.0E0
            IDI(I)=1
            ICOUNT=ICOUNT+1
          ENDIF
        ENDDO
        IF(ICOUNT.EQ.3)THEN
          IADJUST=1
        ELSEIF(ICOUNT.EQ.2)THEN
          IADJUST=1
          DO I=1,3
            K1=K21(1,I)
            K2=K21(2,I)
            K3=K21(3,I)
            IF(IDI(K1).NE.-1 .AND. IDI(K1).NE.1)THEN
              DI(K1)=1.0E0-DI(K2)-DI(K3)
              IDI(K1)=1
              IF(DABS(DI(K1)).LE.DN_SAFE)IDI(K1)=-1
              GOTO 350
            ENDIF
          ENDDO
        ELSEIF(ICOUNT.EQ.1)THEN
          IADJUST=1
          DO I=1,3
            K1=K21(1,I)
            K2=K21(2,I)
            K3=K21(3,I)
            IF(IDI(K1).EQ.-1 .OR. IDI(K1).EQ.1)THEN
              DIFF=1.0E0-DI(K1)-DI(K2)-DI(K3)
              DI(K2)=DI(K2)+DIFF
              IF(DABS(DI(K2)).LE.DN_SAFE)THEN
                DI(K2)=0.0E0
                IDI(K2)=-1
                DI(K3)=1.0E0-DI(K1)-DI(K2)
                IDI(K3)=1
                IF(DABS(DI(K3)).LE.DN_SAFE)IDI(K3)=-1
                GOTO 350
              ELSEIF(DABS(DI(K2)-1.0E0).LE.DN_SAFE)THEN
                DI(K2)=1.0E0
                IDI(K2)=1
                DI(K3)=0.0E0
                IDI(K3)=-1
                GOTO 350
              ENDIF
            ENDIF
          ENDDO
        ENDIF               
C
  350   CONTINUE
        DN(1)=0.5E0*(1.0E0-XI(1))*DI(1)
        DN(2)=0.5E0*(1.0E0-XI(1))*DI(2)
        DN(3)=0.5E0*(1.0E0-XI(1))*DI(3)
        DN(4)=0.5E0*(1.0E0+XI(1))*DI(1)
        DN(5)=0.5E0*(1.0E0+XI(1))*DI(2)
        DN(6)=0.5E0*(1.0E0+XI(1))*DI(3)
C
C
C ===== FOR THE CASE OF A 3-D HEXAHEDRAL ELEMENT
C
      ELSEIF(NEQ.EQ.3 .AND. NODE.EQ.8)THEN
        IXI(1)=0
        IF(DABS(XI(1)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=1.0E0
          IXI(1)=1
        ELSEIF(DABS(XI(1)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(1)=-1.0E0
          IXI(1)=-1
        ENDIF
        IXI(2)=0
        IF(DABS(XI(2)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(2)=1.0E0
          IXI(2)=1
        ELSEIF(DABS(XI(2)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(2)=-1.0E0
          IXI(2)=-1
        ENDIF
        IXI(3)=0
        IF(DABS(XI(3)-1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(3)=1.0E0
          IXI(3)=1
        ELSEIF(DABS(XI(3)+1.0E0).LE.DN_SAFE2)THEN
          IADJUST=1
          XI(3)=-1.0E0
          IXI(3)=-1
        ENDIF
        DN(1)=0.125E0*(1.0E0-XI(1))*(1.0E0-XI(2))*(1.0E0-XI(3))
        DN(2)=0.125E0*(1.0E0+XI(1))*(1.0E0-XI(2))*(1.0E0-XI(3))
        DN(3)=0.125E0*(1.0E0+XI(1))*(1.0E0+XI(2))*(1.0E0-XI(3))
        DN(4)=0.125E0*(1.0E0-XI(1))*(1.0E0+XI(2))*(1.0E0-XI(3))
        DN(5)=0.125E0*(1.0E0-XI(1))*(1.0E0-XI(2))*(1.0E0+XI(3))
        DN(6)=0.125E0*(1.0E0+XI(1))*(1.0E0-XI(2))*(1.0E0+XI(3))
        DN(7)=0.125E0*(1.0E0+XI(1))*(1.0E0+XI(2))*(1.0E0+XI(3))
        DN(8)=0.125E0*(1.0E0-XI(1))*(1.0E0+XI(2))*(1.0E0+XI(3))
      ENDIF
C
  999 CONTINUE
      RETURN
      END
C      
C
C
      SUBROUTINE XSI_2 
     I    (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O     XSI,ETA)
C
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE >
C COMPUTE THE LOCAL COORDINATES BASED ON THE GIVEN CARTESIAN 
C COORDINATES
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
C
      DATA NITER /100/
C
      A1=(XW(1,1)+XW(1,2)+XW(1,3)+XW(1,4))*0.25D0
      A2=(-XW(1,1)+XW(1,2)+XW(1,3)-XW(1,4))*0.25D0
      A3=(-XW(1,1)+XW(1,4)-XW(1,2)+XW(1,3))*0.25D0
      A4=(XW(1,1)-XW(1,4)-XW(1,2)+XW(1,3))*0.25D0
      B1=(XW(2,1)+XW(2,2)+XW(2,3)+XW(2,4))*0.25D0
      B2=(-XW(2,1)+XW(2,2)+XW(2,3)-XW(2,4))*0.25D0
      B3=(-XW(2,1)+XW(2,4)-XW(2,2)+XW(2,3))*0.25D0
      B4=(XW(2,1)-XW(2,4)-XW(2,2)+XW(2,3))*0.25D0
C
C MAKE THE INITIAL GUESS OF XSI AND ETA
C
      XSIO=0.0E0
      ETAO=0.0E0
C
C START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C === COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=XTEMP(1)-A1-A2*XSIO-A3*ETAO-A4*XSIO*ETAO
        F2=XTEMP(2)-B1-B2*XSIO-B3*ETAO-B4*XSIO*ETAO
        Z11=-A2-A4*ETAO
        Z12=-A3-A4*XSIO
        Z21=-B2-B4*ETAO
        Z22=-B3-B4*XSIO
C
C === SOLVE FOR DELXSI AND DELETA
C
        DJAC=Z11*Z22 - Z21*Z12
        DELXSI=(F1*Z22-F2*Z12)/DJAC
        DELETA=(-F1*Z21+F2*Z11)/DJAC
C
C === COMPUTE FOR NEW XSI AND ETA
C
        XSI=XSIO-DELXSI
        ETA=ETAO-DELETA
C
C === TEST CONVERGENCE
C
        D1=DELXSI*0.5E0*DL
        D2=DELETA*0.5E0*DL
        IF(DABS(D1).LE.ATOL .AND. DABS(D2).LE.ATOL)THEN
          RETURN
        ENDIF
C
C === UPDATE XSIO AND ETAO
C
        XSIO=XSI
        ETAO=ETA
  590 CONTINUE
C
C === CONVERGENCE FAILS, STOP EXECUTION
C
      WRITE(*,*)'WARNING IN XSI_2'
      WRITE(*,2000)ITER,NITER,D1,D2,DL,ATOL
 2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING XSI, ETA:',/1X,
     >       'ITER =',I3,',  NITER =',I3,',  D1 =',E15.6/1X,
     >       'D2 =',E15.6,',  DL=',E15.6,',   ATOL =',E15.6)
      WRITE(*,*)'XSI,ETA=',XSI,ETA
      WRITE(*,*)'XTEMP(1:2)=',XTEMP(1),XTEMP(2)
      WRITE(*,*)'XW(1:2,1:4) ='
      WRITE(*,*)((XW(I,J),I=1,2),J=1,4)
      STOP
C
C ===== CONVERGENT SOLUTION FOR XSI AND ETA HAS BEEN ACHIEVED.
C
  999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE XSI_3
     I   (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O    XSI,ETA,ZTA)
C
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE >
C COMPUTE THE LOCAL COORDINATES BASED ON THE GIVEN CARTESIAN 
C COORDINATES
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
C
      DATA NITER/100/
C
C
C MAKE AN INITIAL GUESS OF XSI, ETA, AND ZTA
C
      XSIO=0.0E0
      ETAO=0.0E0
      ZTAO=0.0E0
C
      A1=(XW(1,1)+XW(1,2)+XW(1,3)+XW(1,4)+XW(1,5)+XW(1,6)+
     >    XW(1,7)+XW(1,8))*0.125E0
      A2=(-XW(1,1)+XW(1,2)+XW(1,3)-XW(1,4)-XW(1,5)+XW(1,6)+
     >    XW(1,7)-XW(1,8))*0.125E0
      A3=(-XW(1,1)+XW(1,4)-XW(1,2)+XW(1,3)-XW(1,5)-XW(1,6)+
     >    XW(1,7)+XW(1,8))*0.125E0
      A4=(XW(1,1)-XW(1,4)-XW(1,2)+XW(1,3)+XW(1,5)-XW(1,6)+
     >    XW(1,7)-XW(1,8))*0.125E0
      A5=(-XW(1,1)-XW(1,2)-XW(1,3)-XW(1,4)+XW(1,5)+XW(1,6)+
     >    XW(1,7)+XW(1,8))*0.125E0
      A6=(XW(1,1)+XW(1,2)-XW(1,3)-XW(1,4)-XW(1,5)-XW(1,6)+
     >    XW(1,7)+XW(1,8))*0.125E0
      A7=(XW(1,1)-XW(1,2)-XW(1,3)+XW(1,4)-XW(1,5)+XW(1,6)+
     >    XW(1,7)-XW(1,8))*0.125E0
      A8=(-XW(1,1)+XW(1,2)-XW(1,3)+XW(1,4)+XW(1,5)-XW(1,6)+
     >    XW(1,7)-XW(1,8))*0.125E0
      B1=(XW(2,1)+XW(2,2)+XW(2,3)+XW(2,4)+XW(2,5)+XW(2,6)+
     >    XW(2,7)+XW(2,8))*0.125E0
      B2=(-XW(2,1)+XW(2,2)+XW(2,3)-XW(2,4)-XW(2,5)+XW(2,6)+
     >    XW(2,7)-XW(2,8))*0.125E0
      B3=(-XW(2,1)+XW(2,4)-XW(2,2)+XW(2,3)-XW(2,5)-XW(2,6)+
     >    XW(2,7)+XW(2,8))*0.125E0
      B4=(XW(2,1)-XW(2,4)-XW(2,2)+XW(2,3)+XW(2,5)-XW(2,6)+
     >    XW(2,7)-XW(2,8))*0.125E0
      B5=(-XW(2,1)-XW(2,2)-XW(2,3)-XW(2,4)+XW(2,5)+XW(2,6)+
     >    XW(2,7)+XW(2,8))*0.125E0
      B6=(XW(2,1)+XW(2,2)-XW(2,3)-XW(2,4)-XW(2,5)-XW(2,6)+
     >    XW(2,7)+XW(2,8))*0.125E0
      B7=(XW(2,1)-XW(2,2)-XW(2,3)+XW(2,4)-XW(2,5)+XW(2,6)+
     >    XW(2,7)-XW(2,8))*0.125E0
      B8=(-XW(2,1)+XW(2,2)-XW(2,3)+XW(2,4)+XW(2,5)-XW(2,6)+
     >    XW(2,7)-XW(2,8))*0.125E0
      C1=(XW(3,1)+XW(3,2)+XW(3,3)+XW(3,4)+XW(3,5)+XW(3,6)+
     >    XW(3,7)+XW(3,8))*0.125E0
      C2=(-XW(3,1)+XW(3,2)+XW(3,3)-XW(3,4)-XW(3,5)+XW(3,6)+
     >    XW(3,7)-XW(3,8))*0.125E0
      C3=(-XW(3,1)+XW(3,4)-XW(3,2)+XW(3,3)-XW(3,5)-XW(3,6)+
     >    XW(3,7)+XW(3,8))*0.125E0
      C4=(XW(3,1)-XW(3,4)-XW(3,2)+XW(3,3)+XW(3,5)-XW(3,6)+
     >    XW(3,7)-XW(3,8))*0.125E0
      C5=(-XW(3,1)-XW(3,2)-XW(3,3)-XW(3,4)+XW(3,5)+XW(3,6)+
     >    XW(3,7)+XW(3,8))*0.125E0
      C6=(XW(3,1)+XW(3,2)-XW(3,3)-XW(3,4)-XW(3,5)-XW(3,6)+
     >    XW(3,7)+XW(3,8))*0.125E0
      C7=(XW(3,1)-XW(3,2)-XW(3,3)+XW(3,4)-XW(3,5)+XW(3,6)+
     >    XW(3,7)-XW(3,8))*0.125E0
      C8=(-XW(3,1)+XW(3,2)-XW(3,3)+XW(3,4)+XW(3,5)-XW(3,6)+
     >    XW(3,7)-XW(3,8))*0.125E0
C
C
C START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C === COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=XTEMP(1)-A1-A2*XSIO-A3*ETAO-A4*XSIO*ETAO-A5*ZTAO-
     >              A6*ETAO*ZTAO-A7*ZTAO*XSIO-A8*XSIO*ETAO*ZTAO
        F2=XTEMP(2)-B1-B2*XSIO-B3*ETAO-B4*XSIO*ETAO-B5*ZTAO-
     >              B6*ETAO*ZTAO-B7*ZTAO*XSIO-B8*XSIO*ETAO*ZTAO
        F3=XTEMP(3)-C1-C2*XSIO-C3*ETAO-C4*XSIO*ETAO-C5*ZTAO-
     >              C6*ETAO*ZTAO-C7*ZTAO*XSIO-C8*XSIO*ETAO*ZTAO
        Z11=-A2-A4*ETAO-A7*ZTAO-A8*ETAO*ZTAO
        Z12=-A3-A4*XSIO-A6*ZTAO-A8*ZTAO*XSIO
        Z13=-A5-A6*ETAO-A7*XSIO-A8*XSIO*ETAO
        Z21=-B2-B4*ETAO-B7*ZTAO-B8*ETAO*ZTAO
        Z22=-B3-B4*XSIO-B6*ZTAO-B8*ZTAO*XSIO
        Z23=-B5-B6*ETAO-B7*XSIO-B8*XSIO*ETAO
        Z31=-C2-C4*ETAO-C7*ZTAO-C8*ETAO*ZTAO
        Z32=-C3-C4*XSIO-C6*ZTAO-C8*ZTAO*XSIO
        Z33=-C5-C6*ETAO-C7*XSIO-C8*XSIO*ETAO
C
C === SOLVE FOR DELXSI, DELETA, AND DELZT
C
        DJAC = Z11*(Z22*Z33-Z23*Z32)-Z21*(Z12*Z33-Z32*Z13)+
     1         Z31*(Z12*Z23-Z22*Z13)
        DELXSI=(F1*(Z22*Z33-Z32*Z23)-F2*(Z12*Z33-Z32*Z13)+
     1         F3*(Z12*Z23-Z22*Z13))/DJAC
        DELETA=(-F1*(Z21*Z33-Z31*Z23)+F2*(Z11*Z33-Z31*Z13)-
     1         F3*(Z11*Z23-Z21*Z13))/DJAC
        DELZTA=(F1*(Z21*Z32-Z31*Z22)-F2*(Z11*Z32-Z31*Z12)+
     1         F3*(Z11*Z22-Z21*Z12))/DJAC
C
C === COMPUTE FOR NEW XSI, ETA, AND ZTA
C
        XSI=XSIO-DELXSI
        ETA=ETAO-DELETA
        ZTA=ZTAO-DELZTA
C
C === EXAMINE CONVERGENCE
C
        D1=DELXSI*0.5E0*DL
        D2=DELETA*0.5E0*DL
        D3=DELZTA*0.5E0*DL
        IF(DABS(D1).LE.ATOL .AND. DABS(D2).LE.ATOL .AND. 
     >     DABS(D3).LE.ATOL)THEN
          RETURN
        ENDIF
C
C === UPDATE XSIO, ETAO, AND ZTAO
C
        XSIO=XSI
        ETAO=ETA
        ZTAO=ZTA
C
  590 CONTINUE
C
C
C CONVERGENCE FAILS, STOP EXECUTION
C
      WRITE(*,*)'WARNING IN XSI_3'
      WRITE(*,2000)ITER,NITER,D1,D2,D3,DL,ATOL
 2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING DL1, DL2, XSI:',/1X,
     >       'ITER =',I3,',  NITER =',I3,/1X,
     >       'D1 =',E15.6/1X,'  D2 =',E15.6,',  DL3 =',E15.6,1X,
     >       'DL=',E15.6,',  ATOL =',E15.6)
      WRITE(*,*)'XSI,ETA,ZTA=',XSI,ETA,ZTA
      WRITE(*,*)'XTEMP(1:3)=',XTEMP(1),XTEMP(2),XTEMP(3)
      WRITE(*,*)'XW(1:3,1:8) ='
      DO J=1,8
        WRITE(*,*)(XW(I,J),I=1,3)
      ENDDO
      STOP
C
C CONVERGENT SOLUTION FOR XSI, ETA, AND ZTA HAS BEEN ACHIEVED.
C
  999 CONTINUE
      RETURN
      END
C
C
C 
      SUBROUTINE XSI_3P
     I   (MAXEQ,MAXND,XTEMP,XW,DL,ATOL,
     O    XSI,DL1,DL2,DL3)
C
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE >
C COMPUTE THE LOCAL COORDINATES BASED ON THE GIVEN CARTESIAN 
C COORDINATES
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XTEMP(MAXEQ),XW(MAXEQ,MAXND)
C
      DATA NITER/100/
C
C
C MAKE AN INITIAL GUESS OF XSIO, L1O, L2O
C
      XSIO=0.0E0
      DL1O=0.33333333E0
      DL2O=0.33333333E0
C
C COMPUTE COOEFFICIENTS TO BE USED IN NEWTON-RALPHSON METHOD
C
      A1=XW(1,4)+XW(1,1)
      A2=XW(1,5)+XW(1,2)
      A3=XW(1,6)+XW(1,3)
      A4=XW(1,4)-XW(1,1)
      A5=XW(1,5)-XW(1,2)
      A6=XW(1,6)-XW(1,3)
      B1=XW(2,4)+XW(2,1)
      B2=XW(2,5)+XW(2,2)
      B3=XW(2,6)+XW(2,3)
      B4=XW(2,4)-XW(2,1)
      B5=XW(2,5)-XW(2,2)
      B6=XW(2,6)-XW(2,3)
      C1=XW(3,4)+XW(3,1)
      C2=XW(3,5)+XW(3,2)
      C3=XW(3,6)+XW(3,3)
      C4=XW(3,4)-XW(3,1)
      C5=XW(3,5)-XW(3,2)
      C6=XW(3,6)-XW(3,3)
C
C
C START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C === COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=0.5E0*((A1-A3)*DL1O+(A2-A3)*DL2O)+
     >     0.5E0*((A4-A6)*DL1O+(A5-A6)*DL2O)*XSIO+
     >     0.5E0*A6*XSIO-(XTEMP(1)-0.5E0*A3)
        F2=0.5E0*((B1-B3)*DL1O+(B2-B3)*DL2O)+
     >     0.5E0*((B4-B6)*DL1O+(B5-B6)*DL2O)*XSIO+
     >     0.5E0*B6*XSIO-(XTEMP(2)-0.5E0*B3)
        F3=0.5E0*((C1-C3)*DL1O+(C2-C3)*DL2O)+
     >     0.5E0*((C4-C6)*DL1O+(C5-C6)*DL2O)*XSIO+
     >     0.5E0*C6*XSIO-(XTEMP(3)-0.5E0*C3)
C
        Z11=0.5E0*(A1-A3)+0.5E0*(A4-A6)*XSIO
        Z12=0.5E0*(A2-A3)+0.5E0*(A5-A6)*XSIO
        Z13=0.5E0*((A4-A6)*DL1O+(A5-A6)*DL2O)+0.5E0*A6
        Z21=0.5E0*(B1-B3)+0.5E0*(B4-B6)*XSIO
        Z22=0.5E0*(B2-B3)+0.5E0*(B5-B6)*XSIO
        Z23=0.5E0*((B4-B6)*DL1O+(B5-B6)*DL2O)+0.5E0*B6
        Z31=0.5E0*(C1-C3)+0.5E0*(C4-C6)*XSIO
        Z32=0.5E0*(C2-C3)+0.5E0*(C5-C6)*XSIO
        Z33=0.5E0*((C4-C6)*DL1O+(C5-C6)*DL2O)+0.5E0*C6
C
C === SOLVE FOR DELXSI, DELETA, AND DELZT
C
        DJAC=Z11*(Z22*Z33-Z23*Z32)-Z21*(Z12*Z33-Z32*Z13) +
     >       Z31*(Z12*Z23-Z22*Z13)
        DELDL1=(F1*(Z22*Z33-Z32*Z23)-F2*(Z12*Z33-Z32*Z13)+
     >          F3*(Z12*Z23-Z22*Z13))/DJAC
        DELDL2=(-F1*(Z21*Z33-Z31*Z23)+F2*(Z11*Z33-Z31*Z13)-
     >          F3*(Z11*Z23-Z21*Z13))/DJAC
        DELXSI=(F1*(Z21*Z32-Z31*Z22)-F2*(Z11*Z32-Z31*Z12)+
     >          F3*(Z11*Z22-Z21*Z12))/DJAC
C
C === COMPUTE FOR NEW XSI, ETA, AND ZTA
C
        DL1=DL1O-DELDL1
        DL2=DL2O-DELDL2
        XSI=XSIO-DELXSI
C
C === EXAMINE CONVERGENCE
C
        D1=DELDL1*DL
        D2=DELDL2*DL
        D3=DELXSI*0.5E0*DL
        IF(DABS(D1).LE.ATOL .AND. DABS(D2).LE.ATOL .AND. 
     >     DABS(D3).LE.ATOL)THEN
          GOTO 999
        ENDIF
C
C === UPDATE DL1O, DL2O, AND XSIO
C
        DL1O=DL1
        DL2O=DL2
        XSIO=XSI
C
  590 CONTINUE
C
C
C CONVERGENCE FAILS, STOP EXECUTION
C
      WRITE(*,*)'WARNING IN XSI_3P'
      WRITE(*,2000)ITER,NITER,D1,D2,D3,DL,ATOL
 2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING DL1, DL2, XSI:',/1X,
     >       'ITER =',I3,',  NITER =',I3,/1X,
     >       'D1 =',E15.6/1X,'  D2 =',E15.6,',  DL3 =',E15.6,1X,
     >       'DL=',E15.6,',  ATOL =',E15.6)
      WRITE(*,*)'DL1,DL2,XSI=',DL1,DL2,XSI
      WRITE(*,*)'XTEMP(1:3)=',XTEMP(1),XTEMP(2),XTEMP(3)
      WRITE(*,*)'XW(1:3,1:6) ='
      DO J=1,6
        WRITE(*,*)(XW(I,J),I=1,3)
      ENDDO
      STOP
C
C
C CONVERGENT SOLUTION FOR XSI, ETA, AND ZTA HAS BEEN ACHIEVED.
C
  999 CONTINUE
      DL3=1.0E0-DL1-DL2
C
      RETURN
      END
C
C
C
      SUBROUTINE CL123
     I    (MAXEQ,MAXND,NEQ,NODE,XW,
     O     DL)
C
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C   COMPUTE THE CHARACTERISTIC LENGTH OF A GIVEN ELEMENT
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION KVB3(4,6,3),KVB2(3,2,2)
      DIMENSION XX(4),YY(4),ZZ(4)
      DATA KVB3 /5,1,2,4, 6,5,2,4, 8,5,6,4, 7,2,3,4, 8,6,7,4, 6,7,4,2,
     >           4,1,2,3, 5,4,2,3, 6,4,5,3, 0,0,0,0, 0,0,0,0, 0,0,0,0,
     >           4,1,2,3, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0/
      DATA KVB2 /1,2,3, 1,3,4,
     >           1,2,3, 0,0,0/
C
C
C ===== FOR THE CASE OF A 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
        DL=DABS(XW(1,2)-XW(1,1))
C
C
C ===== FOR THE CASE OF A 2-D ELEMENT
C
      ELSEIF(NEQ.EQ.2)THEN
C
C ... FOR A TRIANGULAR ELEMENT
C
        IF(NODE.EQ.3)THEN
          ID=2
          NA=1
C
C ... FOR A QUADRILATER ELEMENT
C
        ELSE
          ID=1
          NA=2
        ENDIF
C
        DJAC=0.0D0
        DO IA=1,NA
          DO I=1,3
            II=KVB2(I,IA,ID)
            XX(I)=XW(1,II)
            YY(I)=XW(2,II)
          ENDDO
          X12=XX(1)-XX(2)
          X23=XX(2)-XX(3)
          X31=XX(3)-XX(1)  
          Y12=YY(1)-YY(2)
          Y23=YY(2)-YY(3)
          Y31=YY(3)-YY(1)
          DJAC=DJAC+XX(1)*Y23+XX(2)*Y31+XX(3)*Y12
        ENDDO
        DJAC=DABS(DJAC)/DBLE(NA)
        DL=DSQRT(DJAC)
C
C
C ===== FOR THE CASE OF A 3-D ELEMENT
C
      ELSEIF(NEQ.EQ.3)THEN
C
C ... FOR A TETRAHEDRAL ELEMENT
C
        IF(NODE.EQ.4)THEN
          ID=3
          NV=1
C
C ... FOR A TRIANGULAR PRISM ELEMENT
C
        ELSEIF(NODE.EQ.6)THEN
          ID=2
          NV=3
C
C ... FOR A HEXAHEDRAL ELEMENT
C
        ELSE
          ID=1
          NV=6
        ENDIF
C
        DJAC=0.0D0
        DO IV=1,NV
          DO I=1,4
            II=KVB3(I,IV,ID)
            XX(I)=XW(1,II)
            YY(I)=XW(2,II)
            ZZ(I)=XW(3,II)
          ENDDO
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
            DJAC=DJAC+
     >           (-1.0D0)**(KK+1)*(XX(K1)*YY(K2)*ZZ(K3)+
     >            YY(K1)*ZZ(K2)*XX(K3)+ZZ(K1)*XX(K2)*YY(K3)-
     >            XX(K3)*YY(K2)*ZZ(K1)-YY(K3)*ZZ(K2)*XX(K1)-
     >            ZZ(K3)*XX(K2)*YY(K1))
          ENDDO
        ENDDO
        DJAC=DABS(DJAC)/DBLE(NV)
        DL=DJAC**(1.0E0/3.0E0)
      ENDIF
C
  999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE PHI_COMP
     I    (MAXND,NODE,NEQ, 
     I     DN_SAFE,  
     I     DN_S,DN,XI,XI_S,DI,DI_S, IXI_S,IDI_S,
     O     IDSDT,I1,I2,I3,PHI)
C
C 07/14/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C LOOP OVER ALL POSSIBLE SCENARIOS TO (1) DETERMINE PHI IF THE PT PASSES 
C THROUGH THE ELEMENT AND ENDS AT A LOCATION OUTSIDE OF THE ELEMENT
C AND (2) LOCATE I1, I2, I3 TO CONTINUE PARTICLE TRACKING WHEN THE PT
C DOES NOT PASS THROUGH THIS ELEMENT
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DN_S(MAXND),DN(MAXND)
      DIMENSION XI(3),XI_S(3),DI(3),DI_S(3)
      DIMENSION IXI_S(3),IDI_S(4)
C
      DN_SAFE2=2.0E0*DN_SAFE
      IF(NEQ.EQ.2 .AND. NODE.NE.3)GOTO 200
      IF(NEQ.EQ.3 .AND. NODE.NE.4)GOTO 300
C
C === COMPUTE PHI
C
C (1) WHEN 1-D LINE ELEMENT, 2-D TRIANGULAR ELEMENT, 
C     OR 3-D TETRAHEDRAL ELEMENT IS CONSIDERED
C
  100 CONTINUE
C
      DO 150 I=1,NODE
C
C ... WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C
        IF(DN(I).LT.0.0E0)THEN
          D1=DN_S(I)
          D2=DN(I)
          D12=D1-D2
          CALL PHI123
     I        (NEQ,NODE, IXI_S,IDI_S, 
     I         D1,D2,D12, DN_SAFE,  
     O         IDSDT,I1,I2,I3,PHI)
          IF(IDSDT.EQ.0)RETURN
        ENDIF
  150 CONTINUE
      RETURN
C
C (2) WHEN 2-D QUADRILATERAL ELEMENT IS CONSIDERED
C
  200 CONTINUE
C
      DO 250 I=1,2
C
C ... WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C
        IF(XI(I).LT.-1.0E0)THEN
          D1=XI_S(I)+1.0E0
          D2=XI(I)+1.0E0
          D12=D1-D2
          CALL PHI123
     I        (NEQ,NODE, IXI_S,IDI_S, 
     I         D1,D2,D12, DN_SAFE2,  
     O         IDSDT,I1,I2,I3,PHI)
          IF(IDSDT.EQ.0)RETURN
        ENDIF
        IF(XI(I).GT.1.0E0)THEN
          D1=1.0E0-XI_S(I)
          D2=1.0E0-XI(I)
          D12=D1-D2
          CALL PHI123
     I        (NEQ,NODE, IXI_S,IDI_S, 
     I         D1,D2,D12, DN_SAFE2, 
     O         IDSDT,I1,I2,I3,PHI)
          IF(IDSDT.EQ.0)RETURN
        ENDIF
  250 CONTINUE
      RETURN
C
C (3) WHEN 3-D TRIANGULAR PRISM OR HEXAHETRAL ELEMENT
C     IS CONSIDERED
C
  300 CONTINUE
C
C (3.1) FOR THE CASE OF TRIANGULAR PRISM ELEMENT
C
      IF(NODE.EQ.6)THEN
C
C ... WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C
        IF(XI(1).LT.-1.0E0)THEN
          D1=XI_S(1)+1.0E0
          D2=XI(1)+1.0E0
          D12=D1-D2
          CALL PHI123
     I        (NEQ,NODE, IXI_S,IDI_S, 
     I         D1,D2,D12, DN_SAFE2, 
     O         IDSDT,I1,I2,I3,PHI)
          IF(IDSDT.EQ.0)RETURN
        ENDIF
        IF(XI(1).GT.1.0E0)THEN
          D1=1.0E0-XI_S(1)
          D2=1.0E0-XI(1)
          D12=D1-D2
          CALL PHI123
     I        (NEQ,NODE, IXI_S,IDI_S, 
     I         D1,D2,D12, DN_SAFE2, 
     O         IDSDT,I1,I2,I3,PHI)
          IF(IDSDT.EQ.0)RETURN
        ENDIF
C
        DO 350 I=1,3
C
C ... WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C
          IF(DI(I).LT.0.0E0)THEN
            D1=DI_S(I)
            D2=DI(I)
            D12=D1-D2
            CALL PHI123
     I          (NEQ,NODE, IXI_S,IDI_S, 
     I           D1,D2,D12, DN_SAFE, 
     O           IDSDT,I1,I2,I3,PHI)
            IF(IDSDT.EQ.0)RETURN
          ENDIF
  350   CONTINUE
C
C (3.2) FOR THE CASE OF HEXAHEDRAL ELEMENT
C
      ELSE
        DO 370 I=1,3
C
C ... WHEN THE ENDING LOCATION IS OUTSIDE OF ELEMENT M
C
          IF(XI(I).LT.-1.0E0)THEN
            D1=XI_S(I)+1.0E0
            D2=XI(I)+1.0E0
            D12=D1-D2
            CALL PHI123
     I          (NEQ,NODE, IXI_S,IDI_S, 
     I           D1,D2,D12, DN_SAFE2, 
     O           IDSDT,I1,I2,I3,PHI)
            IF(IDSDT.EQ.0)RETURN
          ENDIF
          IF(XI(I).GT.1.0E0)THEN
            D1=1.0E0-XI_S(I)
            D2=1.0E0-XI(I)
            D12=D1-D2
            CALL PHI123
     I          (NEQ,NODE, IXI_S,IDI_S, 
     I           D1,D2,D12, DN_SAFE2, 
     O           IDSDT,I1,I2,I3,PHI)
            IF(IDSDT.EQ.0)RETURN
          ENDIF
  370   CONTINUE
      ENDIF
C
  999 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE PHI123
     I    (NEQ,NODE, IXI_S,IDI_S,
     I     D1,D2,D12, D_SAFE, 
     O     IDSDT,I1,I2,I3,PHI)
C 
C 07/14/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C COMPUTE PHI BASED ON THE GIVEN D1, D2, AND D12
C LOCATE I1, I2, I3 IF THE PARTICLE DOES NOT PASS THROUGH THE ELEMENT
C < NOTE >
C D1 = THE DISTANCE OF THE STARTING LOCATION FROM A SPECIFIED ELEMENT 
C      BOUNDARY
C D2 = THE DISTANCE OF THE ENDING LOCATION FROM THE SAME SPECIFIED 
C      ELEMENT BOUNDARY
C D12 = THE DISTANCE BETWEEN THE STARTING AND THE ENDING LOCATIONS 
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IXI_S(3),IDI_S(4)
C
C (1) FOR THE CASE WHEN PT IS NOT THROUGH THIS ELEMENT
C
      IF(DABS(D1).LT.D_SAFE)THEN
C
C IF THE TRACKING IS LEAVING THE ELEMENT
C ==> SET IDSDT TO 0
C     IDENTIFY I1, I2, I3 FOR FINDING AN ADJACENT ELEMENT TO CONTINUE PT
C
C < NOTE > IF THE TRACKING IS ALONG AN ELEMENT BOUNDARY 
C          ==> NO ACTION IS NEEDED
C
        IF(DABS(D2).GT.D_SAFE)THEN
          IDSDT=0
          CALL EB_CHECK
     I        (NEQ,NODE, IDI_S,IXI_S, 
     O         I1,I2,I3)

          RETURN
        ENDIF
C
C (2) FOR THE CASE WHEN PT IS THROUGH THIS ELEMENT
C
      ELSE
C
C IF THE ENDING LOCATION CAN BE CONSIDERED ON THE ELEMENT BOUNDARY
C ==> NO NEED TO REDUCE TIMESTEP
C
        IF(DABS(D2).LE.D_SAFE)RETURN
C
C ... IF THE ENDING LOCATION IS TRULY OUTSIDE OF THE ELEMENT
C     ==> COMPUTE AN INTERPOLATION FACTOR
C
        PHI=DMAX1(PHI,D12/D1)
      ENDIF
C
  999 CONTINUE
      RETURN
      END
C
C 
C
      SUBROUTINE EB_CHECK
     I    (NEQ,NODE, IDI_S,IXI_S, 
     O     I1,I2,I3)
C 
C 07/14/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C PREPARE I1, I2, I3 FOR THE SUCCESSIVE PT WHEN THE TRACKED NODE IS
C ON THE ELEMENT BOUNDARY
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IXI_S(3),IDI_S(4)
C
C ===== INITIALIZATION
C
      I1=0
      I2=0
      I3=0
C
C ===== IDENTIFY I1, I2, I3
C
C (1) 1-D LINE ELEMENT
C
      IF(NEQ.EQ.1)THEN
        IF(IDI_S(1).EQ.1)THEN
          I1=1
          RETURN
        ELSEIF(IDI_S(2).EQ.1)THEN
          I1=2
          RETURN
        ENDIF
C
C (2) 2-D ELEMENT
C
      ELSEIF(NEQ.EQ.2)THEN
C
C ... FOR TRIANGULAR ELEMENT
C
        IF(NODE.EQ.3)THEN
          DO I=1,NODE
            IF(IDI_S(I).EQ.1)THEN
              I1=I
              RETURN
            ENDIF
          ENDDO
C
          IF(IDI_S(1).EQ.-1)THEN
            I1=2
            I2=3
            RETURN
          ELSEIF(IDI_S(2).EQ.-1)THEN
            I1=1
            I2=3
            RETURN
          ELSEIF(IDI_S(3).EQ.-1)THEN
            I1=1
            I2=2
            RETURN
          ENDIF
C
C ... FOR QUADRILATERAL ELEMENT
C
        ELSEIF(NODE.EQ.4)THEN
          IF(IXI_S(1).EQ.-1)THEN
            IF(IXI_S(2).EQ.-1)THEN
              I1=1
              RETURN
            ELSEIF(IXI_S(2).EQ.1)THEN
              I1=4
              RETURN
            ELSE
              I1=1
              I2=4
              RETURN
            ENDIF
          ELSEIF(IXI_S(1).EQ.1)THEN
            IF(IXI_S(2).EQ.-1)THEN
              I1=2
              RETURN
            ELSEIF(IXI_S(2).EQ.1)THEN
              I1=3
              RETURN
            ELSE
              I1=2
              I2=3
              RETURN
            ENDIF
          ELSE
            IF(IXI_S(2).EQ.-1)THEN
              I1=1
              I2=2
              RETURN
            ELSEIF(IXI_S(2).EQ.1)THEN
              I1=3
              I2=4
              RETURN
            ENDIF
          ENDIF  
        ENDIF 
C
C (3) 3-D ELEMENT
C  
      ELSEIF(NEQ.EQ.3)THEN
C
C ... FOR TETRAHEDRAL ELEMENT
C
        IF(NODE.EQ.4)THEN  
          DO I=1,NODE
            IF(IDI_S(I).EQ.1)THEN
              I1=I
              RETURN
            ENDIF
          ENDDO
C
          IF(IDI_S(1).EQ.-1)THEN
            IF(IDI_S(2).EQ.-1)THEN
              I1=3
              I2=4
              RETURN
            ELSEIF(IDI_S(3).EQ.-1)THEN
              I1=2
              I2=4
              RETURN
            ELSEIF(IDI_S(4).EQ.-1)THEN
              I1=2
              I2=3
              RETURN
            ELSE
              I1=2
              I2=3
              I3=4
              RETURN
            ENDIF
          ELSEIF(IDI_S(2).EQ.-1)THEN
            IF(IDI_S(3).EQ.-1)THEN
              I1=1
              I2=4
              RETURN
            ELSEIF(IDI_S(4).EQ.-1)THEN
              I1=1
              I2=3
              RETURN
            ELSE
              I1=1
              I2=3
              I3=4
              RETURN
            ENDIF            
          ELSEIF(IDI_S(3).EQ.-1)THEN
            IF(IDI_S(4).EQ.-1)THEN
              I1=1
              I2=2
              RETURN
            ELSE
              I1=1
              I2=2
              I3=4
              RETURN
            ENDIF
          ELSEIF(IDI_S(4).EQ.-1)THEN
            I1=1
            I2=2
            I3=3
            RETURN
          ENDIF              
C
C ... FOR TRIANGULAR PRISM ELEMENT
C
        ELSEIF(NODE.EQ.6)THEN 
          IF(IXI_S(1).EQ.-1)THEN
            DO I=1,3
              IF(IDI_S(I).EQ.1)THEN
                I1=I
                RETURN
              ENDIF
            ENDDO
            IF(IDI_S(1).EQ.-1)THEN
              I1=2
              I2=3
            ELSEIF(IDI_S(2).EQ.-1)THEN
              I1=1
              I2=3
            ELSEIF(IDI_S(3).EQ.-1)THEN
              I1=1
              I2=2
            ELSE
              I1=1
              I2=2
              I3=3
            ENDIF
            RETURN
          ELSEIF(IXI_S(1).EQ.1)THEN
            DO I=1,3
              IF(IDI_S(I).EQ.1)THEN
                I1=I+3
                RETURN
              ENDIF
            ENDDO
            IF(IDI_S(1).EQ.-1)THEN
              I1=5
              I2=6
            ELSEIF(IDI_S(2).EQ.-1)THEN
              I1=4
              I2=6
            ELSEIF(IDI_S(3).EQ.-1)THEN
              I1=4
              I2=5
            ELSE
              I1=4
              I2=5
              I3=6
            ENDIF
            RETURN
          ELSE
            DO I=1,3
              IF(IDI_S(I).EQ.1)THEN
                I1=I
                I2=I+3
                RETURN
              ENDIF
            ENDDO
            IF(IDI_S(1).EQ.-1)THEN
              I1=2
              I2=3
              I3=5
              RETURN
            ELSEIF(IDI_S(2).EQ.-1)THEN
              I1=1
              I2=3
              I3=6
              RETURN
            ELSEIF(IDI_S(3).EQ.-1)THEN
              I1=1
              I2=2
              I3=4
              RETURN
            ENDIF
          ENDIF
C
C ... FOR HEXAHEDRAL ELEMENT
C
        ELSEIF(NODE.EQ.8)THEN  
          IF(IXI_S(1).EQ.-1)THEN
            IF(IXI_S(2).EQ.-1)THEN
              IF(IXI_S(3).EQ.-1)THEN
                I1=1
                RETURN
              ELSEIF(IXI_S(3).EQ.1)THEN
                I1=5
                RETURN
              ELSE
                I1=1
                I2=5
              ENDIF
            ELSEIF(IXI_S(2).EQ.1)THEN
              IF(IXI_S(3).EQ.-1)THEN
                I1=4
                RETURN
              ELSEIF(IXI_S(3).EQ.1)THEN
                I1=8
                RETURN
              ELSE
                I1=4
                I2=8
              ENDIF
            ELSEIF(IXI_S(3).EQ.-1)THEN
              I1=1
              I2=4
              RETURN
            ELSEIF(IXI_S(3).EQ.1)THEN
              I1=5
              I2=8
              RETURN
            ELSE
              I1=1
              I2=4
              I3=5
            ENDIF
          ELSEIF(IXI_S(1).EQ.1)THEN
            IF(IXI_S(2).EQ.-1)THEN
              IF(IXI_S(3).EQ.-1)THEN
                I1=2
                RETURN
              ELSEIF(IXI_S(3).EQ.1)THEN
                I1=6
                RETURN
              ELSE
                I1=2
                I2=6
              ENDIF
            ELSEIF(IXI_S(2).EQ.1)THEN
              IF(IXI_S(3).EQ.-1)THEN
                I1=3
                RETURN
              ELSEIF(IXI_S(3).EQ.1)THEN
                I1=7
                RETURN
              ELSE
                I1=3
                I2=7
              ENDIF
            ELSEIF(IXI_S(3).EQ.-1)THEN
              I1=2
              I2=3
              RETURN
            ELSEIF(IXI_S(3).EQ.1)THEN
              I1=6
              I2=7
              RETURN
            ELSE
              I1=2
              I2=3
              I3=6
            ENDIF
          ELSEIF(IXI_S(2).EQ.-1)THEN
            IF(IXI_S(3).EQ.-1)THEN
              I1=1
              I2=2
              RETURN
            ELSEIF(IXI_S(3).EQ.1)THEN
              I1=5
              I2=6
              RETURN
            ELSE
              I1=1
              I2=2
              I3=5
            ENDIF
          ELSEIF(IXI_S(2).EQ.1)THEN
            IF(IXI_S(3).EQ.-1)THEN
              I1=3
              I2=4
              RETURN
            ELSEIF(IXI_S(3).EQ.1)THEN
              I1=7
              I2=8
              RETURN
            ELSE
              I1=3
              I2=4
              I3=7
            ENDIF
          ELSEIF(IXI_S(3).EQ.-1)THEN
            I1=1
            I2=2
            I3=3
            RETURN
          ELSEIF(IXI_S(3).EQ.1)THEN
            I1=5
            I2=6
            I3=7
            RETURN
          ENDIF  
        ENDIF
C
      ENDIF
C
  999 CONTINUE
      RETURN
      END 
C
C 
C
      SUBROUTINE DT_ETRACK
     I    (MAXND,MAXEQ,NODE,NEQ, 
     I     DL,VT1W,VT2W,
     O     DT_INIT)
C 
C 09/09/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C ESTIMATE  ELEMENT NODAL VELOCITY FOR PT WITHIN THE ELEMENT
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
C 
      VABS=0.0D0
      NCOUNT=NODE*NEQ*2
      DO I=1,NODE
        DO J=1,NEQ   
          VABS=VABS+DABS(VT1W(J,I)) 
          VABS=VABS+DABS(VT2W(J,I))
        ENDDO
      ENDDO
      VABS=VABS/DBLE(NCOUNT)
      IF(VABS.LT.1.0E-10)VABS=1.0E-10
      DT_INIT=DL/VABS 
C
  999 CONTINUE
      RETURN
      END

C ======================================================================
CMWF ADD ROUTINES FOR RT0 WITH NUMERICAL TRACKING
C ======================================================================
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
      SUBROUTINE EL_VEL_PREP
     I    (MAXND,MAXEQ,NNDE,
     I     NNP,NEL,NODE,NEQ, M,IDVE,DIR,IPROJ,
     I     XG,IE,IB,VTL2G,VT1E,VT2E,
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
C TRY TO PROJECT ON BOUNDARY OR NOT
      INTEGER IPROJ
C NODAL COORDINATES
      DOUBLE PRECISION XG(MAXEQ,NNP)
C ELEMENT -- NODES CONNECTIVITY
      INTEGER IE(NNDE,NEL)
C BOUNDARY NODE INFORMATION: 1 --> ON BOUNDARY, 0 --> INTERIOR
      INTEGER IB(*)
C LOCAL TO GLOBAL DOF MAPPING
C      INTEGER VTL2G(NEQ,NNDE,NEL)
      INTEGER VTL2G(*)
C DEGREES OF FREEDOM AT TIME LEVELS 1 AND 2
      DOUBLE PRECISION VT1E(*),VT2E(*)
C WORK ARRAY FOR NODES AND VELOCITIES AT TIME LEVEL 1 AND 2
      DOUBLE PRECISION XW(MAXEQ,MAXND),VT1W(MAXEQ,MAXND),
     &     VT2W(MAXEQ,MAXND)
C LOCAL VARIABLES
      INTEGER I,J,K,IVDOF,IEM,ID,NS,KS,NCOUNT,N1,N2,N3,N4,K1,K2,K3,K4
      DOUBLE PRECISION DVOL,DRT0FACT
C WORK ARRAYS FOR LOCAL PROJECTION
      INTEGER KGB3(4,6,3),KGB2(2,4,2)
C
      DATA KGB3 /1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >           1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >           4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
      DATA KGB2 /1,2, 2,3, 3,4, 4,1,
     >           1,2, 2,3, 3,1, 0,0/
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
C            WRITE(6,*)'EL_VEL_PREP IDVE.EQ.2 M= ',M,' J= ',J,
C     &           ' IEM= ',IEM,' XW(',K,',',J,')= ',XW(K,J),
C     &           ' IVDOF= ',IVDOF,' VT1W(',K,',',J,')= ',VT1W(K,J),
C     &           ' VT2W(',K,',',J,')= ',VT2W(K,J),' IPROJ= ',IPROJ
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
C
CMWF ADDED CHECK TO MAKE SURE NOT 1D
      IF(IPROJ.EQ.0.OR.NEQ.EQ.1)RETURN
C
C === CONDUCT VELOCITY PROJECTION ONTO CLOSED BOUNDARY WHEN NECESSARY
C
C 2-D PROJECTION
C
      IF(NEQ.EQ.2)THEN
C
C FOR A TRIANGULAR ELEMENT
        IF(NODE.EQ.3)THEN
          ID=2
          NS=3
C FOR A QUADRILATERAL ELEMENT
        ELSE
          ID=1
          NS=4
        ENDIF    
C
        DO 100 KS=1,NS
          K1=KGB2(1,KS,ID)
          K2=KGB2(2,KS,ID)
          N1=IE(K1,M)
          N2=IE(K2,M)
          NCOUNT=0
CMWF CONVENTION 1 --> BOUNDARY, PEARCE'S IS -1 --> BOUNDARY
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1)THEN
            CALL V_PROJTN23
     I          (MAXND,MAXEQ,
     I           NEQ,K1,K2,K3,K4,
     M           XW,VT1W,VT2W)
          ENDIF
  100   CONTINUE
C
C 3-D PROJECTION
C
      ELSE
C
C FOR A TETRAHEDRAL ELEMENT
        IF(NODE.EQ.4)THEN
          ID=3
          NS=4
C FOR A TRIANGULAR PRISM ELEMENT
        ELSEIF(NODE.EQ.6)THEN
          ID=2
          NS=5
C FOR A HEXAHEDRAL ELEMENT
        ELSE
          ID=1
          NS=6
        ENDIF
C 
        DO 200 KS=1,NS
          K1=KGB3(1,KS,ID)
          K2=KGB3(2,KS,ID)
          K3=KGB3(3,KS,ID)
          K4=KGB3(4,KS,ID)
          N1=IE(K1,M)
          N2=IE(K2,M)
          N3=IE(K3,M)
CMWF CONVENTION 1 --> BOUNDARY, PEARCE'S IS -1 --> BOUNDARY
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1)THEN
            IF(K4.NE.0)THEN
              N4=IE(K4,M)
              IF(IB(N4).NE.-1)GOTO 200
            ENDIF
            CALL V_PROJTN23
     I          (MAXND,MAXEQ,
     I           NEQ,K1,K2,K3,K4,
     M           XW,VT1W,VT2W)
          ENDIF
  200   CONTINUE
C
      ENDIF      
C

      RETURN
      END
C
C
C
      SUBROUTINE V_PROJTN23
     I    (MAXND,MAXEQ,
     I     NEQ,K1,K2,K3,K4,
     M     XW,VT1W,VT2W)
C 
C 06/22/2010 (HPC)
C ======================================================================
C < PURPOSE > 
C CONDUCT VELOCITY PROJECT ONTO A 2-D ELEMENT EDGE OR A 3-D ELEMENT FACE
C IF IT IS ASSOCIATED WITH CLOSED BOUNDARY
C ======================================================================
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XW(MAXEQ,MAXND)
      DIMENSION VT1W(MAXEQ,MAXND),VT2W(MAXEQ,MAXND)
C
C ===== FOR PROJECTION ONTO A 2-D ELEMENET EDGE COMPOSED OF 
C       ELEMENT NODES K1 AND K2
C
      IF(NEQ.EQ.2)THEN
        DX=XW(1,K2)-XW(1,K1)
        DY=XW(2,K2)-XW(2,K1)
        DL=DSQRT(DX*DX+DY*DY)
        DL2=DL*DL
C NODE K1:
        PK1T1=(VT1W(1,K1)*DX+VT1W(2,K1)*DY)/DL2
        VT1W(1,K1)=PK1T1*DX
        VT1W(2,K1)=PK1T1*DY  
        PK1T2=(VT2W(1,K1)*DX+VT2W(2,K1)*DY)/DL2
        VT2W(1,K1)=PK1T2*DX
        VT2W(2,K1)=PK1T2*DY
C NODE K2:   
        PK2T1=(VT1W(1,K2)*DX+VT1W(2,K2)*DY)/DL2
        VT1W(1,K2)=PK2T1*DX
        VT1W(2,K2)=PK2T1*DY  
        PK2T2=(VT2W(1,K2)*DX+VT2W(2,K2)*DY)/DL2
        VT2W(1,K2)=PK2T2*DX
        VT2W(2,K2)=PK2T2*DY     
C
C ===== FOR PROJECTION UNTO A 3-D ELEMENET FACE COMPOSED OF 
C       ELEMENT NODES K1, K2, K3
C
      ELSEIF(NEQ.EQ.3)THEN
        DX12=XW(1,K2)-XW(1,K1)
        DY12=XW(2,K2)-XW(2,K1)
        DZ12=XW(3,K2)-XW(3,K1)
        DX13=XW(1,K3)-XW(1,K1)
        DY13=XW(2,K3)-XW(2,K1)
        DZ13=XW(3,K3)-XW(3,K1)
        DNX=DY12*DZ13-DY13*DZ12
        DNY=DZ12*DX13-DZ13*DX12
        DNZ=DX12*DY13-DX13*DY12
        DLN=DSQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
        DLN2=DLN*DLN
C NODE K1:
        PK1T1N=(VT1W(1,K1)*DNX+VT1W(2,K1)*DNY+VT1W(3,K1)*DNZ)/DLN2
        VT1W(1,K1)=VT1W(1,K1)-PK1T1N*DNX   
        VT1W(2,K1)=VT1W(2,K1)-PK1T1N*DNY  
        VT1W(3,K1)=VT1W(3,K1)-PK1T1N*DNZ  
        PK1T2N=(VT2W(1,K1)*DNX+VT2W(2,K1)*DNY+VT2W(3,K1)*DNZ)/DLN2
        VT2W(1,K1)=VT2W(1,K1)-PK1T2N*DNX   
        VT2W(2,K1)=VT2W(2,K1)-PK1T2N*DNY  
        VT2W(3,K1)=VT2W(3,K1)-PK1T2N*DNZ
C NODE K2:
        PK2T1N=(VT1W(1,K2)*DNX+VT1W(2,K2)*DNY+VT1W(3,K2)*DNZ)/DLN2
        VT1W(1,K2)=VT1W(1,K2)-PK2T1N*DNX   
        VT1W(2,K2)=VT1W(2,K2)-PK2T1N*DNY  
        VT1W(3,K2)=VT1W(3,K2)-PK2T1N*DNZ  
        PK2T2N=(VT2W(1,K2)*DNX+VT2W(2,K2)*DNY+VT2W(3,K2)*DNZ)/DLN2
        VT2W(1,K2)=VT2W(1,K2)-PK2T2N*DNX   
        VT2W(2,K2)=VT2W(2,K2)-PK2T2N*DNY  
        VT2W(3,K2)=VT2W(3,K2)-PK2T2N*DNZ    
C NODE K3:
        PK3T1N=(VT1W(1,K3)*DNX+VT1W(2,K3)*DNY+VT1W(3,K3)*DNZ)/DLN2
        VT1W(1,K3)=VT1W(1,K3)-PK3T1N*DNX   
        VT1W(2,K3)=VT1W(2,K3)-PK3T1N*DNY  
        VT1W(3,K3)=VT1W(3,K3)-PK3T1N*DNZ  
        PK3T2N=(VT2W(1,K3)*DNX+VT2W(2,K3)*DNY+VT2W(3,K3)*DNZ)/DLN2
        VT2W(1,K3)=VT2W(1,K3)-PK3T2N*DNX   
        VT2W(2,K3)=VT2W(2,K3)-PK3T2N*DNY  
        VT2W(3,K3)=VT2W(3,K3)-PK3T2N*DNZ  
C NODE K4:
        IF(K4.NE.0)THEN
          PK4T1N=(VT1W(1,K4)*DNX+VT1W(2,K4)*DNY+VT1W(3,K4)*DNZ)/DLN2
          VT1W(1,K4)=VT1W(1,K4)-PK4T1N*DNX   
          VT1W(2,K4)=VT1W(2,K4)-PK4T1N*DNY  
          VT1W(3,K4)=VT1W(3,K4)-PK4T1N*DNZ  
          PK4T2N=(VT2W(1,K4)*DNX+VT2W(2,K4)*DNY+VT2W(3,K4)*DNZ)/DLN2
          VT2W(1,K4)=VT2W(1,K4)-PK4T2N*DNX   
          VT2W(2,K4)=VT2W(2,K4)-PK4T2N*DNY  
          VT2W(3,K4)=VT2W(3,K4)-PK4T2N*DNZ 
        ENDIF
      ENDIF
C
  999 CONTINUE
      RETURN
      END

C ======================================================================
C Start working on analytical tracking in separate routine, then merge
C ======================================================================

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
      IF (IVERBOSE.GT.0) THEN
         WRITE(LU_O,*)'ENTERING PT123A'
         WRITE(LU_O,*)'NNP= ',NNP,' NEL= ',NEL,' NNDE= ',NNDE,' NEQ= ',
     &     NEQ,' NPT= ',NPT,' IBF= ',IBF,' LU_O= ',LU_O, ' ATOL= ',
     &     ATOL,' RTOL= ',RTOL, ' SF= ',SF, 'DN_SAFE= ',DN_SAFE
      ENDIF
      DO IPT=1,NPT
         IF(IVERBOSE.GT.1) THEN
CMWF DEBUG
            WRITE(LU_O,*)'T_I(',IPT,')= ',T_I(IPT),
     &           ' TPT(',IPT,')= ',TPT(IPT)
         ENDIF
         DO I=1,NEQ
            XS(I)=X_I(I + MAXEQ*(IPT-1))
            XPT(I + MAXEQ*(IPT-1))=XS(I)
            IF(IVERBOSE.GT.1) THEN
CMWF DEBUG
               WRITE(LU_O,*)'X_I(',I,',',IPT,')= ',X_I(I+MAXEQ*(IPT-1)),
     &              ' XPT(',I,',',IPT,')= ',XPT(I + MAXEQ*(IPT-1))
            ENDIF
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
          CALL EL_VEL_PREPA(MAXND,MAXEQ,NNDE,
     &         NNP,NEL,NODE,NEQ,M,IDVE,DIR,
     &         XG,IE,VTL2G,VT1E,VT2E,XW,VT1W,VT2W)

C
C CONDUCT TRACKING WITHIN ELEMENT M
C
          CALL ELTRAK123ANEW
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
        IF(IVERBOSE.GT.1) THEN
C MWF DEBUG
           WRITE(LU_O,*)' B4 ELTRACK ELEMENT LOOP T= ',T,' TPT= ',
     &          TPT(IPT),' IPT= ',IPT,' M= ',M
           WRITE(LU_O,*)' XS= ',(XS(I),I=1,NEQ)
        ENDIF
        CALL EL_VEL_PREPA(MAXND,MAXEQ,NNDE,
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
        CALL ELTRAK123ANEW
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
        IF(IVERBOSE.GT.1) THEN
C MWF DEBUG
           WRITE(LU_O,*)' AFTER ELTRACK ELEMENT LOOP T= ',T,
     &          ' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT,' M= ',M
           WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
           DO K=1,NEQ
              WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',XPT(K+(IPT-1)*MAXEQ)
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
                   CALL EL_VEL_PREPA(MAXND,MAXEQ,NNDE,
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
                  CALL ELTRAK123ANEW
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
                  IF(IVERBOSE.GT.2) THEN
C MWF DEBUG
                     WRITE(LU_O,*)'AFTER I3.NE.0 ELTRACK ELEMENT LOOP',
     &                    ' T= ',T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),
     &                    ' IPT= ',IPT,' M= ',M
                     WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
                     DO K=1,NEQ
                        WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &                       XPT(K+(IPT-1)*MAXEQ)
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
                 CALL EL_VEL_PREPA(MAXND,MAXEQ,NNDE,
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
                CALL ELTRAK123ANEW
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
                IF(IVERBOSE.GT.2) THEN
C MWF DEBUG
                   WRITE(LU_O,*)' AFTER I2.NE.0 ELTRACK ELEMENT LOOP',
     &               ' T= ',T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),
     &               ' IPT= ',IPT,' M= ',M
                   WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
                   DO K=1,NEQ
                      WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &                     XPT(K+(IPT-1)*MAXEQ)
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
          CALL EL_VEL_PREPA(MAXND,MAXEQ,NNDE,
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
          CALL ELTRAK123ANEW
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
          IF(IVERBOSE.GT.2) THEN
C MWF DEBUG
             WRITE(LU_O,*)' AFTER I1.NE.0 ELTRACK ELEMENT LOOP T= '
     &            ,T,' IDSDT= ',IDSDT,' TPT= ',TPT(IPT),' IPT= ',IPT
     &         ,' M= ',M
             WRITE(LU_O,*)' I1= ',I1,' I2= ',I2,' I3= ',I3
             DO K=1,NEQ
                WRITE(LU_O,*)' XPT(',K,',',IPT,')= ',
     &               XPT(K+(IPT-1)*MAXEQ)
             ENDDO
          ENDIF
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
C RT0, LINEAR DEPENDENCE IN TIME
         CALL STEPDTRT0V1
     I        (MAXEQ,MAXND,NEQ,NODE,
     I        DEQ, DN_SAFE,
     I        ZEROTOL,DIR,T1,T2,
     I        XW,VT1W,VT2W,
     I        FNORMALS,XFBAR,
     M        T,DT,XS,XOUT5,XTEMP,
     M        VTEMP,DN_S,DN,DDN)
      ELSE IF (IVFLAG.EQ.2) THEN
C RT0, BILINEAR IN X-T
         CALL STEPDTRT0V2
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
      CALL INTRP123A
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
            IF(DABS(D1).LT.DN_SAFE)THEN
C            IF(DABS(D1*DL).LT.ZEROTOL)THEN
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
C PYADH HAD               IF(DABS(D2).GT.DN_SAFE)THEN
C                  CALL DN_CHECKA
C     I                 (MAXND,NODE,NEQ, DL,ZEROTOL,
C     I                 DN_S,
C     O                 I1,I2,I3)
               IF(DABS(D2*DL).GT.ZEROTOL)THEN
                  CALL DN_CHECKA
     I                 (MAXND,NODE,NEQ, DL,ZEROTOL,
     I                 DN_S,
     O                 I1,I2,I3)
                  IDSDT=0

                  IF (IDEBUG.GT.0) THEN
C MWF DEBUG
                     WRITE(6,*)' ELTRAK DN_CHECKA LEAVING IDSDT= ',
     &                    IDSDT,' I1= ',I1,' I2= ',I2,' I3= ',I3, 
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
C PYADH HAD               IF(DABS(D2).LE.DN_SAFE) THEN
               IF(DABS(D2).LE.DN_SAFE) THEN
C LATER VERSION WAS
C               IF(DABS(D2*DL).LE.ZEROTOL) THEN
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
C PYADH HAD
C         IF(DABS(DN(I)).LT.DN_SAFE)THEN
C            CALL DN_CHECK
C     I           (MAXND,NODE,NEQ, 1.D0,DN_SAFE,
C     I           DN,
C     O           I1,I2,I3)
C LATER VERSION WAS
C         IF(DABS(DN(I)*DL).LT.ZEROTOL)THEN
C            CALL DN_CHECKA
C     I           (MAXND,NODE,NEQ, DL,ZEROTOL,
C     I           DN,
C     O           I1,I2,I3)
         IF(DABS(DN(I)*DL).LT.ZEROTOL)THEN
            CALL DN_CHECKA
     I           (MAXND,NODE,NEQ, DL,ZEROTOL,
     I           DN,
     O           I1,I2,I3)
            IF (IDEBUG.GT.0) THEN
C MWF DEBUG
               WRITE(6,*)' ELTRAK DN_CHECKA ON BNDY RETUR IDSDT= ',
     &              IDSDT,' I1= ',I1,' I2= ',I2,' I3= ',I3, 
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
            WRITE(6,*)' ELTRAK DN_CHECKA ON BNDY RETURN IDSDT= ',IDSDT, 
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

C
C
C
CMWF ORIGINAL VERSION OF INTRP123 ROUTINE
      SUBROUTINE INTRP123A
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
CMWF PREVIOUS VERSION OF EL_VEL_PREP USED RIGHT NOW IN ELTRAK123A
      SUBROUTINE EL_VEL_PREPA
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
C            WRITE(6,*)'EL_VEL_PREP IDVE.EQ.2 M= ',M,' J= ',J,
C     &           ' IEM= ',IEM,' XW(',K,',',J,')= ',XW(K,J),
C     &           ' IVDOF= ',IVDOF,' VT1W(',K,',',J,')= ',VT1W(K,J),
C     &           ' VT2W(',K,',',J,')= ',VT2W(K,J)
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
C
C 
C
      SUBROUTINE DN_CHECKA
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
C
C
C 
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
      CALL INTRP123A
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
C THIS VERSION USES CUTOFF CONSIDERATIONS FROM HEALY AND RUSSELL EQNS (16)--(19)
      SUBROUTINE DXUPDATE1
     I     (ZEROTOL,DT,V0,VX,VT,
     O     DX)
C
      IMPLICIT NONE
C EPSILON, MACHINE_EPS
      DOUBLE PRECISION ZEROTOL
C \Delta t, v_0, v_x, v_t, v_{xt}
      DOUBLE PRECISION DT,V0,VX,VT,VXT
C x-x_0
      DOUBLE PRECISION DX
C LOCAL VARIABLES
C DELTA CUTOFF FOR v_xt AND S=v_x\Delta t FROM FORMULA 18 AND 19
      DOUBLE PRECISION DELTA,S,SMAG,SE
C     \frac{e^s-1}{v_x}, \frac{e^s-(1+s)}{v_x^2}
      DOUBLE PRECISION TERM1,TERM2
C     \frac{1}{v_x}\left(\frac{1}{2}\Delta t^2 e^s - \frac{e^s-(1+s)}{v^2_x}\right)
      DOUBLE PRECISION TERM3
C     \frac{e^s - (1+s+\frac{1}{2}s^2 + \frac{1}{6}s^3}{v_x^4}
      DOUBLE PRECISION TERM4
C     \frac{e^s - (1+s+\frac{1}{2}s^2)}{v^3_x}
      DOUBLE PRECISION TERM5
C     
      DOUBLE PRECISION ONETHIRD,ONESIXTH,ONEEIGHTH,ONE24

      ONETHIRD = 1.D0/3.D0
      ONESIXTH = 1.D0/6.D0
      ONEEIGHTH= 1.D0/8.D0
      ONE24    = 1.D0/24.D0
      S    = VX*DT
      SMAG = DABS(S)
      SE   = DEXP(S)
C EVALUATE TERMS ON PAGE 6, RH00
      IF (SMAG.LT.(6.D0*ZEROTOL)**ONETHIRD) THEN
         TERM1 = DT*(1.D0 + 0.5*S)
      ELSE
         TERM1 = (SE-1.D0)/VX
      ENDIF
      IF (SMAG.LT.(24.D0*ZEROTOL)**(0.25D0)) THEN
         TERM2 = DT*DT*(0.5 + S*ONESIXTH)
      ELSE
         TERM2 = (SE-(1.D0+S))/(VX*VX)
      ENDIF
      IF (SMAG.LT.(40.D0*ONETHIRD*ZEROTOL)**(0.2D0)) THEN
         TERM3 = (DT**3)*(ONETHIRD + 5.D0*S*ONE24)
      ELSE
         TERM3 = (0.5*DT*DT*SE - (SE - (1.D0+S))/(VX*VX))/VX
      ENDIF
      IF (SMAG.LT.(720.D0*ZEROTOL)**ONESIXTH) THEN
         TERM4 = DT**(4.D0)*(ONE24 + S*ONESIXTH*0.05D0)
      ELSE
         TERM4 = (SE - (1.D0 + S + 0.5*S**2 + ONESIXTH*S**3))/(VX**4)
      ENDIF
      IF (SMAG.LT.(120.D0*ZEROTOL)*(0.2D0)) THEN
         TERM5 = DT**(3.D0)*(ONESIXTH + ONE24*S)
      ELSE
         TERM5 = (SE - (1.D0 + S + 0.5*S*S))/(VX**3)
      ENDIF

      DX = V0*TERM1 + VT*TERM2

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
      DOUBLE PRECISION DX0F,V0F,VTF,VX,DXF
C     DOUBLE PRECISION ALPHA,DALPHA,BETA,DBETA
      DOUBLE PRECISION ZEROTOL
      
      DX0F= XPAR(1)
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)
      
C      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
C      R = DX0F + ALPHA*V0F + BETA*VTF
C      DR= DALPHA*V0F + DBETA*VTF
      CALL DXUPDATE1(ZEROTOL,DT,V0F,VX,VTF,DXF)
      R = DX0F + DXF
      DR= V0F + VX*DXF + VTF*DT

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
      DOUBLE PRECISION VF,DXF,V0F,VTF,VX
C     DOUBLE PRECISION ALPHA,DALPHA,BETA,DBETA
      DOUBLE PRECISION ZEROTOL
            
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)

C      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
C      DXF = ALPHA*V0F + BETA*VTF
      CALL DXUPDATE1(ZEROTOL,DT,V0F,VX,VTF,DXF)
      VF  = V0F + VX*DXF + VTF*DT
      R   = VF
      DR  = VX*VF + VTF

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
      DOUBLE PRECISION X1(3),DX10(3),DXOUT(3),
     &     V1(3),VT(3),VTEMP2(3),RPAR(10)
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
      MAXIT = 200
C 1 - FORCE CONVERGENCE ON RESIDUAL AND DX 
C 0 - ALLOW EITHER
      IPAR(1) = 1
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
      CALL INTRP123A
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
      IF (IPAR(1).EQ.0) THEN
C BE STRICTER IF ALLOWING CONVERGENCE ONLY ON |NEWTON INCR|
         DXINCRTOL= 1.D-3*DXNLTOL
      ELSE
         DXINCRTOL=0.1D0*DXNLTOL
      ENDIF
C MAKE CONSISTENT WITH TESTS WHEN VELOCITY IS POSTITIVE OR NEGATIVE
      VNLTOL  = ZEROTOL
      IF (IPAR(1).EQ.0) THEN
C BE STRICTER IF ALLOWING CONVERGENCE ONLY ON |NEWTON INCR|
         VINCRTOL= 1.D-3*VNLTOL 
      ELSE
         VINCRTOL= 0.1D0*VNLTOL 
      ENDIF
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
C      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
      DO I=1,NEQ
C         X1(I) = XS(I) + ALPHA*VTEMP(I) + BETA*VT(I)
         CALL DXUPDATE1(ZEROTOL,DT,VTEMP(I),VX,VT(I),DX10(I))
         X1(I) = XS(I) + DX10(I)
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
C TODO HAVE TO CHECK THIS REASONING ON DIR AND 
C      DECISION WITH MAGNITUDE OF VTF
            IF (DABS(VTF).GE.1.D-3*ZEROTOL) THEN
               SV0F = DSIGN(1.D0,VTF*DIR)
            ELSE
               SV0F = DIR
            ENDIF
         ELSE
            SV0F = DSIGN(1.D0,V0F)
         ENDIF
         IF (DABS(V1F).LE.1.D-3*ZEROTOL) THEN
C  TODO HAVE TO CHECK THIS REASONING ON DIR AND
C      DECISION WITH MAGNITUDE OF VTF1
            IF (DABS(VTF).GE.1.D-3*ZEROTOL) THEN
               SV1F = DSIGN(1.D0,-VTF*DIR)
            ELSE
               SV1F = DIR
            ENDIF
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
C                  CALL DXCOEFS1(ZEROTOL,DTFI,VX,ALPHA,DALPHA,BETA,DBETA)
C                  DXFI = DXF0 + ALPHA*V0F + BETA*VTF
                  CALL DXUPDATE1(ZEROTOL,DTFI,V0F,VX,VTF,DXFI)
                  DXFI = DXF0 + DXFI
                  WRITE(6,*)'STEPDTRT0V1 CASE A DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFI= ',DTFI,' DTF= ',DTF,' DXFI= ',DXFI
                  WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &                 ' DL= ',DL,'IPAR(3)= ',IPAR(3),
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
C            CALL DXCOEFS1(ZEROTOL,DTFS,VX,ALPHA,DALPHA,BETA,DBETA)
C            DXFS = DXF0 + ALPHA*V0F + BETA*VTF 
            CALL DXUPDATE1(ZEROTOL,DTFS,V0F,VX,VTF,DXFS)
            DXFS = DXF0 + DXFS
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
C               CALL DXCOEFS1(ZEROTOL,DTFS,VX,ALPHA,DALPHA,BETA,DBETA)
C               DXFS = DXF0 + ALPHA*V0F + BETA*VTF 
               CALL DXUPDATE1(ZEROTOL,DTFS,V0F,VX,VTF,DXFS)
               DXFS = DXF0 + DXFS
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
C      CALL DXCOEFS1(ZEROTOL,DT,VX,ALPHA,DALPHA,BETA,DBETA)
      DO I=1,NEQ
C         XOUT(I) = XS(I) + ALPHA*VTEMP(I) + BETA*VT(I)
         CALL DXUPDATE1(ZEROTOL,DT,VTEMP(I),VX,VT(I),DXOUT(I))
         XOUT(I) = XS(I) + DXOUT(I)
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
C BEGIN ROUTINES SPECIFIC TO BILINEAR (X-T) VELOCITY
C
C (CAPITAL LETTERS ARE VECTORS)
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

      SUBROUTINE STEPDTRT0V2
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
C RT0 VELOCITY ON A SIMPLEX WITH FULL BILINIEARITY IN X-T
C 
C RT0 VELOCITY CAN BE REPRESENTED AS (CAPITAL LETTERS ARE VECTORS DIM=NEQ)
C
C  V = V_0 + v_x (X-X_0) + V_t (t-t_0) + v_xt (X-X_0)(t-t_0)
C
C WHERE X_0 IS THE STARTING POSITION WITH VELOCITY V_0=V(X_0)
C THE ANALYTICAL SOLUTION FOR THE POSITION IS GIVEN BY EQUATIONS (12)--(15) IN
C RUSSELL AND HEALY 2000, WITH THE ACTUAL IMPLEMENTATION SPECIFIED IN (16)--(19)
C 
C
C POINT X IS IN SIMPLEX M IFF (X_F-X).N_F <= 0 FOR ALL F IN FACES(M)
C WHERE 
C  X_F IS THE BARYCENTER OF FACE F, WITH UNIT OUTER NORMAL N_F
C 
C THIS APPROACH IS VERY SIMILAR TO THE LINEAR IN X-T CASE, EXCEPT THAT THERE
C CAN BE 2 VELOCITY REVERSALS OVER A TIME STEP IN TWO CASES.
C
C THE DOUBLE REVERSAL CASE OCCURS IFF sign(V_0.N_F) == sign(V_1.N_F) AND
C   sign(V_0.N_F) != sign((X_1-X_0).N_F)
C 
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
C     dx_f10= (X_1 - X_0).N_f
C
C    IF |dx_f0|  <= epsilon, SKIP FACE (ON BOUNDARY OR OUTSIDE ALREADY)
C
C    OTHERWISE EVALUATE CASES BELOW BASED ON v_0f, v_1f, WHERE dir IS THE 
C    DIRECTION GOING IN TIME, sign[v] IS CALCULATED USING v_t IF |v| <= 1.0e-3 epsilon
C    TODO CHECK THIS MODIFICATION OF HEALY AND RUSSELL, BUT I'M USING
C    sign[v_t*dir] IF |v_0f| <= 1.0e-3 epsilon and sign[-v_t*dir] IF |v_1f| <= 1.0e-3 epsilon
C    THEY DO NOT INCLUDE dir IN THE TEST
C
C A. sign[v_f0] == sign[v_f1] == dir and sign[v_0]*sign[dx_f10] == dir (no reversal)
C    if dx_f1  < 0, never intersects face
C    otherwise,
C      solve for timestep  to intersect face, dt_i, using dx_f(dt) = R, bracketed between
C       R(0) = dx_f0 < 0, R(dt_1) = dx_f1 > 0   
C    set dt_out = min(dt_out,dt_i) 
C
C B. sign[v_f0] == sign[v_f1] == -dir and sign[v_0]*sign[dx_f10] == dir (no reversal)
C    never intersects face, skip
C
C C. sign[v_f0] == sign[v_f1] == dir and sign[v_0]*sign[dx_f10] == -dir (two reversals, v_0 towards face)
C    solve for dt_t in [0,dt_1] such that x(dt_t) = (x_0 + x_1)/2, to do this
C       solve for dx_t = x(dt_t)-x_0 = (x_1-x_0)/2 using usual intersection residual routine
C       with dx_f0 replaced by 0.5*dx_01 = (dx_f0-dx_f1)/2, bracket by R(0) and R(dt)
C       Must have dx_t - dx_f0 < 0, or else there is a problem
C 
C    solve for dt_s1 in [0,dt_t] where the first reversal occurs using v_f(dt) = R bracketed
C       by v_0f and v(dt_t),
C    compute dx_s1(dt_s1). 
C
C    If dx_s1 < 0, then no intersection
C    otherwise, 
C      solve for timestep to intersect dt_i, dx_f(dt_i)=0 using dx_f(dt) = R, bracketed between
C      R(0) = dx_0 < 0 and R(dt_s1) = dx_s1 > 0
C      set dt_out = min(dt_out,dt_i)    
C
C D. sign[v_f0] == sign[v_f1] == -dir and sign[v_0]*sign[dx_f10] == -dir (two reversals, v_f0 away from face) 
C    solve for dt_t in [0,dt_1] such that x(dt_t) = (x_0 + x_1)/2, using 
C        R = dx(dt) + (dx_f0-dx_f1)/2, bracketed between R(0) and R(dt)
C        Must have dx_t = (x_1-x_0)/2 > dx_f0 or else there is a problem
C
C    solve for dt_s1 in [0,dt_t] where the first reversal occurs using v_f(dt) = R, bracketed
C       by v_0 = R(0) < 0 and v(dt_t) > 0, (here < 0 means = -dir, > 0 means = dir)
C       Must have dx_s1 = dx(dt_s1) < dx_f0, since v_f0 away from face
C
C    solve for dt_s2 in [dt_t,dt_1] where the second reversal occurs using v_f(dt) = R, bracketed
C       by v_1 = R(dt) < 0 and v(dt_t) > 0, (here < 0 means = -dir, > 0 means = dir)
C
C    if dx_s2 = dx(dt_s2) < 0, then no intersection occurs
C    otherwise, 
C      solve for timestep to intersect face, dt_i, dx_f(dt_i)=0, using dx_f(dt) = R, bracketed between
C      R(0) = dx_s1 < 0 and R(dt_s2) = dx_s2 > 0
C      set dt_out = min(dt_out,dt_1)
C
C E. sign[v_f0] != sign[v_f1], sign[v_f0] = dir
C    solve for dt_s in [0,dt_1], with v_f(dt_s) = 0, using v_f(dt) = R, bracketed by v_0f, v_01
C    compute dx_s(dt_s)
C    if dx_s < 0, never intersects face
C    otherwise,
C      solve for timestep to intersect face, dt_i, dx_f(dt_i)= 0, using dx_f(dt) = R, bracketed between
C       R(0) = dx_f0 < 0, R(dt_s) = dx_s > 0,
C     set dt_out = min(dt_out,dt_i)
C
C F. sign[v_f0] != sign[v_f1], sign[v_f0] = -dir 
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
      EXTERNAL RESVF2,RESDXF2
C LOCAL VARIABLES
      DOUBLE PRECISION X1(3),DX10(3),DXOUT(3),
     &     V1(3),VT(3),VTEMP2(3),RPAR(10)
      INTEGER I,J,II,IDEBUG,NFACE,IPAR(10)
      DOUBLE PRECISION DTF,VTF,VTF1,DXF0,DXF1,DXF10,VX1,VX0,VXT,DXFTT
      DOUBLE PRECISION DTFI,DTFS1,DTFS2,DTFTT,DTOUT,DXFS1,DXFS2,SDXF10
      DOUBLE PRECISION DJAC,DL
      DOUBLE PRECISION V0F,V1F,ALPHA,DALPHA,BETA,DBETA,SV0F,SV1F,DXFI
      DOUBLE PRECISION ETAT0,ETAT1
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
      MAXIT = 200
C FORCE CONVERGENCE ON RESIDUAL AND DX FOR NOW
      IPAR(1) = 1
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V2 BEFORE FACE LOOP T= ',T,' DT= ',DT,
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
      ELSE
         ETAT0 = (T-TV1)/(TV2-TV1)
         ETAT1 = (T+DT-TV1)/(TV2-TV1)
      ENDIF
      CALL INTRP123A
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
      VX0 = 0.D0
      VX1 = 0.D0
      DO J=1,NODE
C ASSUMES RT0 SO DV_XX = DV_YY AND DV_XY=DV_YX=0
         IF (IDEBUG.GT.0) THEN
C MWF DEBUG
            WRITE(6,*)' STEPDTRT0V2 VX0= ',VX0,' VTAVG(1,',J,')= ',
     &           0.5*(VT1W(1,J)+VT2W(1,J)), 'DDN(1,',J,')= ',DDN(1,J)
         ENDIF
C NOW USE AVERAGE OF DV_XX, DV_YY, DV_ZZ
         DO I=1,NEQ
            VX0 = VX0 + DEQ*((1.D0-ETAT0)*VT1W(I,J)*DDN(I,J)+
     &           ETAT0*VT2W(I,J)*DDN(I,J))
            VX1 = VX1 + DEQ*((1.D0-ETAT1)*VT1W(I,J)*DDN(I,J)+
     &           ETAT1*VT2W(I,J)*DDN(I,J))
            VTEMP(I)  = VTEMP(I) + (1.D0-ETAT0)*VT1W(I,J)*DN(J) + 
     &           ETAT0*VT2W(I,J)*DN(J)
            VTEMP2(I) = VTEMP2(I)+ (1.D0-ETAT1)*VT1W(I,J)*DN(J) + 
     &           ETAT1*VT2W(I,J)*DN(J)
         ENDDO
      ENDDO
C MIXED PARTIAL DERIVATIVE
      VXT = (VX1-VX0)/DT
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V2 VX0= ',VX0,' VX1= ',VX1,' VXT= ',VXT,
     &        ' V0= ',(VTEMP(II),II=1,NEQ),
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
      DO I=1,NEQ
         CALL DXUPDATEVXT(ZEROTOL,DT,VTEMP(I),VX0,VT(I),VXT,DX10(I))
         X1(I) = XS(I) + DX10(I)
         V1(I) = VTEMP(I) + VX0*(X1(I)-XS(I)) + VT(I)*DT + 
     &        VXT*DT*(X1(I)-XS(I))
      ENDDO
C
C
CMWF DEBUG
      IF (IDEBUG.GT.0) THEN
         WRITE(6,*)'STEPDTRT0V2 X1= ',(X1(II),II=1,NEQ),
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
         V0F  = 0.D0
         V1F  = 0.D0
         DXF0 = 0.D0
         DXF1 = 0.D0
         DXF10= 0.D0
C INCLUDE V_t.N_f HERE FOR CONVENIENCE, 
         VTF = 0.D0
         VTF1= 0.D0
         DO I=1,NEQ
C NOTE DXF IS OPPOSITE SIGN OF STEADY-STATE EXAMPLE
            V0F  = V0F  + VTEMP(I)*FNORMALS(I,J)
            DXF0 = DXF0 + (XS(I)-XFBAR(I,J))*FNORMALS(I,J)
            V1F  = V1F  + V1(I)*FNORMALS(I,J)
            DXF1 = DXF1 + (X1(I)-XFBAR(I,J))*FNORMALS(I,J)
            VTF  = VTF  + VT(I)*FNORMALS(I,J)
            DXF10= DXF10+ DX10(I)*FNORMALS(I,J)
         ENDDO
C ALSO NEED V_t(x_1,t_1).N_f FOR SIGN CALCULATIONS
         VTF1= VTF + VXT*DXF10
C SIGNS 
         IF (DABS(V0F).LE.1.D-3*ZEROTOL) THEN
C TODO HAVE TO CHECK THIS REASONING ON DIR AND
C      DECISION WITH MAGNITUDE OF VTF
            IF (DABS(VTF).GE.1.D-3*ZEROTOL) THEN
               SV0F = DSIGN(1.D0,VTF*DIR)
            ELSE
               SV0F = DIR
            ENDIF
         ELSE
            SV0F = DSIGN(1.D0,V0F)
         ENDIF
         IF (DABS(V1F).LE.1.D-3*ZEROTOL) THEN
C  TODO HAVE TO CHECK THIS REASONING ON DIR AND
C      DECISION WITH MAGNITUDE OF VTF1
  
            IF (DABS(VTF1).GE.1.D-3*ZEROTOL) THEN
               SV1F = DSIGN(1.D0,-VTF1*DIR)
            ELSE
               SV1F = DIR
            ENDIF
         ELSE
            SV1F = DSIGN(1.D0,V1F)
         ENDIF
C NEED NOW FOR DETERMINING IF HAVE DOUBLE REVERSE
C TODO DETERMINE WHAT'S THE RIGHT THING TO DO IF DXF10 IS SMALL
         IF (DABS(DXF10).GE.1.D-3*ZEROTOL) THEN
            SDXF10 = DSIGN(1.D0,DXF10)
         ELSE
C THIS SHOULD PREVENT DOUBLE REVERSE CASE, OF COURSE NOT SURE IF THAT'S WHAT
C    WE ALWAYS WANT            
            SDXF10 = SV0F*DIR
         ENDIF
         IF (IDEBUG.GT.0) THEN
CMWF DEBUG
            WRITE(6,*)'STEPDTRT0V2 FACE= ',J,' DIR= ',DIR,
     &           ' VTF= ',VTF,' V0F= ',V0F,' SV0F= ',SV0F,
     &           ' V1F= ',V1F,' SV1F= ',SV1F
            WRITE(6,*)'DXF0= ',DXF0,' DXF1= ',DXF1,' DXF10= ',DXF10,
     &           ' SDFX10= ',SDXF10
         ENDIF
C TODO, SWITCH ZEROTOL -> BNDTOL?
C 
         IF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.GT.0.D0.AND.
     &        SV0F*SDXF10*DIR.GT.0.D0) THEN
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
               RPAR(4) = VX0
               RPAR(5) = ZEROTOL
               RPAR(6) = VXT
               DTL = 0.D0
               DTH = DT
               DTFI = RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               DTF  = DTFI
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
C UPDATE POSITION RELATIVE TO FACE
                  CALL DXUPDATEVXT(ZEROTOL,DTFI,V0F,VX0,VTF,VXT,DXFI)
                  DXFI = DXF0 + DXFI
                  WRITE(6,*)'STEPDTRT0V2 CASE A DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFI= ',DTFI,' DTF= ',DTF,' DXFI= ',DXFI
                  WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &                 ' DL= ',DL,'IPAR(3)= ',IPAR(3),
     &                 ' IPAR(2)= ',IPAR(2)
               ENDIF
            ENDIF
         ELSEIF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.LT.0.D0.AND.
     &           SV0F*SDXF10*DIR.GT.0.D0) THEN
C     CASE B
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE B NO INTERSECTION'
            ENDIF
            INTERSECT = .FALSE.
         ELSEIF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.GT.0.D0.AND.
     &           SV0F*SDXF10*DIR.LT.0.D0) THEN
C     CASE C, DOUBLE REVERSE, V0F TOWARDS FACE
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE C DXF10= ',DXF10
            ENDIF
C SOLVE FOR \tilde{dt}, TIME TO REACH MIDPOINT BETWEEN x_1 and x_0 
            DXFTT = DXF10*0.5D0
            RPAR(1) =-DXFTT
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX0
            RPAR(5) = ZEROTOL
            RPAR(6) = VXT
            DTL     = 0.D0
            DTH     = DT
            DTFTT   = RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &           MAXIT,RPAR,IPAR)
C SHOULD HAVE DX(DTFTT) = DXF10/2
            IF (IDEBUG.GT.0) THEN
C COULD DOUBLE CHECK THAT DXF0 + DXFTT = (DXF0+DXF1)/2
               CALL DXUPDATEVXT(ZEROTOL,DTFTT,V0F,VX0,VTF,VXT,DXFTT)
               WRITE(6,*)'STEPDTRT0V2 CASE C DTL= ',DTL,' DTH= ',
     &              DTH,' DTFTT= ',DTFTT,' DXFTT= ',DXFTT,
     &              ' 0.5*DXF10= ',0.5D0*DXF10,' DXF0+DXFTT= ',
     &              DXF0+DXFTT,' 0.5*(DXF0+DXF1)= ',0.5D0*(DXF0+DXF1)
               WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &              ' DL= ',DL,'IPAR(3)= ',IPAR(3),
     &              ' IPAR(2)= ',IPAR(2)
               IF (DABS(0.5D0*DXF10-DXFTT).GT.DXNLTOL) THEN
                  CALL EXIT(1)
               ENDIF
               IF (DABS(DXFTT+DXF0-0.5D0*(DXF0+DXF1)).GT.DXNLTOL) THEN
                  CALL EXIT(1)
               ENDIF
            ENDIF
            IF (DXFTT.GT.DXF0+ZEROTOL) THEN
               WRITE(6,*)'STEPDTRT0V2 CASE C PROBLEM DXFTT= ',DXFTT,
     &              ' SHOULD BE < DXF0= ',DXF0
               CALL EXIT(1)
            ENDIF
C SOLVE FOR FIRST REVERSAL TIME IN [0,\tilde{dt}] = [0,DTFTT]
            RPAR(1) = DXF0
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX0
            RPAR(5) = ZEROTOL
            RPAR(6) = VXT
            DTL  = 0.D0
            DTH  = DTFTT
            DTFS1= RTSAFE(RESVF2,DTL,DTH,VINCRTOL,VNLTOL,MAXIT,
     &           RPAR,IPAR)
C UPDATE POSITION RELATIVE TO FACE
            CALL DXUPDATEVXT(ZEROTOL,DTFS1,V0F,VX0,VTF,VXT,DXFS1)
            DXFS1 = DXF0 + DXFS1
            IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
               WRITE(6,*)'STEPDTRT0V2 CASE C DTL= ',DTL,' DTH= ',
     &              DTH,' DTFS1= ',DTFS1,' DXFS1= ',DXFS1
            ENDIF
            IF (DXFS1.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
C SOLVE FOR INTERSECTION TIME, dt_i in [0,dt_s1]
               RPAR(1) = DXF0
               RPAR(2) = V0F
               RPAR(3) = VTF
               RPAR(4) = VX0
               RPAR(5) = ZEROTOL
               RPAR(6) = VXT
               DTL     = 0.D0
               DTH     = DTFS1
               DTFI = RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
C UPDATE POSITION RELATIVE TO FACE
                  CALL DXUPDATEVXT(ZEROTOL,DTFI,V0F,VX0,VTF,VXT,DXFI)
                  DXFI = DXF0 + DXFI
                  WRITE(6,*)'STEPDTRT0V2 CASE C DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFI= ',DTFI,' DTF= ',DTF,' DXFI= ',DXFI
                  WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &                 ' DL= ',DL,'IPAR(3)= ',IPAR(3),
     &                 ' IPAR(2)= ',IPAR(2)
               ENDIF
               DTF = DTFI
            ENDIF
         ELSEIF (SV0F*SV1F.GT.0.D0.AND.SV0F*DIR.LT.0.D0.AND.
     &           SV0F*SDXF10*DIR.LT.0.D0) THEN
C     CASE D, DOUBLE REVERSE, V0F AWAY FROM FACE
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE D DXF10= ',DXF10
            ENDIF
C SOLVE FOR \tilde{dt}, TIME TO REACH MIDPOINT BETWEEN x_1 and x_0 
            DXFTT = DXF10*0.5D0
            RPAR(1) = -DXFTT
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX0
            RPAR(5) = ZEROTOL
            RPAR(6) = VXT
            DTL     = 0.D0
            DTH     = DT
            DTFTT   = RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &           MAXIT,RPAR,IPAR)
C SHOULD HAVE DX(DTFTT) = DXF10/2
            IF (IDEBUG.GT.0) THEN
C COULD DOUBLE CHECK THAT DXF0 + DXFTT = (DXF0+DXF1)/2
               CALL DXUPDATEVXT(ZEROTOL,DTFTT,V0F,VX0,VTF,VXT,DXFTT)
               WRITE(6,*)'STEPDTRT0V2 CASE D DTL= ',DTL,' DTH= ',
     &              DTH,' DTFTT= ',DTFTT,' DXFTT= ',DXFTT,' DXFBAR= ',
     &              0.5D0*DXF10
               WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &              ' DL= ',DL,'IPAR(3)= ',IPAR(3),
     &              ' IPAR(2)= ',IPAR(2)
               IF (DABS(0.5D0*DXF10-DXFTT).GT.DXNLTOL) THEN
                  CALL EXIT(1)
               ENDIF
               IF (DABS(DXFTT+DXF0-0.5D0*(DXF0+DXF1)).GT.DXNLTOL) THEN
                  CALL EXIT(1)
               ENDIF
            ENDIF
            IF (DXFTT.LT.DXF0+ZEROTOL) THEN
               WRITE(6,*)'STEPDTRT0V2 CASE D PROBLEM DXFTT= ',DXFTT,
     &              ' SHOULD BE > DXF0= ',DXF0
               CALL EXIT(1)
            ENDIF
C SOLVE FOR FIRST REVERSAL TIME IN [0,\tilde{dt}] = [0,DTFTT]
            RPAR(1) = DXF0
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX0
            RPAR(5) = ZEROTOL
            RPAR(6) = VXT
            DTL  = 0.D0
            DTH  = DTFTT
            DTFS1= RTSAFE(RESVF2,DTL,DTH,VINCRTOL,VNLTOL,MAXIT,
     &           RPAR,IPAR)
C UPDATE POSITION RELATIVE TO FACE
            CALL DXUPDATEVXT(ZEROTOL,DTFS1,V0F,VX0,VTF,VXT,DXFS1)
            DXFS1 = DXF0 + DXFS1
            IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
               WRITE(6,*)'STEPDTRT0V2 CASE D DTL= ',DTL,' DTH= ',
     &              DTH,' DTFS1= ',DTFS1,' DXFS1= ',DXFS1
            ENDIF
C SOLVE FOR SECOND REVERSAL TIME IN [\tilde{dt},dt] 
            DTL  = DTFTT
            DTH  = DT
            DTFS2= RTSAFE(RESVF2,DTL,DTH,VINCRTOL,VNLTOL,MAXIT,
     &           RPAR,IPAR)
C UPDATE POSITION RELATIVE TO FACE
            CALL DXUPDATEVXT(ZEROTOL,DTFS2,V0F,VX0,VTF,VXT,DXFS2)
            DXFS2 = DXF0 + DXFS2 
            IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
               WRITE(6,*)'STEPDTRT0V2 CASE D DTL= ',DTL,' DTH= ',
     &              DTH,' DTFS2= ',DTFS2,' DXFS2= ',DXFS2
            ENDIF

            IF (DXFS2.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
C SOLVE FOR INTERSECTION TIME, dt_i in [dt_s1,dt_s2]
               RPAR(1) = DXF0
               RPAR(2) = V0F
               RPAR(3) = VTF
               RPAR(4) = VX0
               RPAR(5) = ZEROTOL
               RPAR(6) = VXT
               DTL     = DTFS1
               DTH     = DTFS2
               DTFI = RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
C UPDATE POSITION RELATIVE TO FACE
                  CALL DXUPDATEVXT(ZEROTOL,DTFI,V0F,VX0,VTF,VXT,DXFI)
                  DXFI = DXF0 + DXFI
                  WRITE(6,*)'STEPDTRT0V2 CASE D DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFI= ',DTFI,' DTF= ',DTF,' DXFI= ',DXFI
                  WRITE(6,*)'DXICRTOL= ',DXINCRTOL,' DXNLTOL= ',DXNLTOL,
     &                 ' DL= ',DL,'IPAR(3)= ',IPAR(3),
     &                 ' IPAR(2)= ',IPAR(2)
               ENDIF
               DTF = DTFI
            ENDIF
         ELSEIF (SV0F*SV1F.LT.0.D0.AND.SV0F*DIR.GT.0.D0) THEN
C     CASE E
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE E DXF1.LT.-BNDTOL= ',DXF1.LT.-BNDTOL
            ENDIF
            RPAR(1) = DXF0
            RPAR(2) = V0F
            RPAR(3) = VTF
            RPAR(4) = VX0
            RPAR(5) = ZEROTOL
            RPAR(6) = VXT
            DTL  = 0.D0
            DTH  = DT
            DTFS1= RTSAFE(RESVF2,DTL,DTH,VINCRTOL,VNLTOL,MAXIT,
     &           RPAR,IPAR)
C UPDATE POSITION RELATIVE TO FACE    
            CALL DXUPDATEVXT(ZEROTOL,DTFS1,V0F,VX0,VTF,VXT,DXFS1)
            DXFS1 = DXF0 + DXFS1
            IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
               WRITE(6,*)'STEPDTRT0V2 CASE E DTL= ',DTL,' DTH= ',
     &              DTH,' DTFS1= ',DTFS1,' DXFS1= ',DXFS1
            ENDIF
            IF (DXFS1.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
               DTL = 0.D0
               DTH = DTFS1
               DTFI= RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &              MAXIT,RPAR,IPAR)
               DTF = DTFI
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
C UPDATE POSITION RELATIVE TO FACE
                  CALL DXUPDATEVXT(ZEROTOL,DTFI,V0F,VX0,VTF,VXT,DXFI)
                  DXFI = DXF0 + DXFI
                  WRITE(6,*)'STEPDTRT0V2 CASE E DTL= ',DTL,
     &                 ' DTH= ',DTH,' DTFI= ',DTFI,' DXFI= ',DXFI,
     &                 ' DTF= ',DTF
               ENDIF
            ENDIF
         ELSEIF (SV0F*SV1F.LT.0.D0.AND.SV0F*DIR.LT.0.D0) THEN
C     CASE F
            IF (IDEBUG.GT.0) THEN
               WRITE(6,*)'CASE F DXF1.LT.-BNDTOL= ',DXF1.LT.-BNDTOL
            ENDIF
            IF (DXF1.LT.-BNDTOL) THEN
               INTERSECT = .FALSE.
            ELSE
               RPAR(1) = DXF0
               RPAR(2) = V0F
               RPAR(3) = VTF
               RPAR(4) = VX0
               RPAR(5) = ZEROTOL
               RPAR(6) = VXT
               DTL  = 0.D0
               DTH  = DT
               DTFS1= RTSAFE(RESVF2,DTL,DTH,VINCRTOL,VNLTOL,
     &              MAXIT,RPAR,IPAR)
C UPDATE POSITION RELATIVE TO FACE
               CALL DXUPDATEVXT(ZEROTOL,DTFS1,V0F,VX0,VTF,VXT,DXFS1)
               DXFS1 = DXF0 + DXFS1
               IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
                  WRITE(6,*)'STEPDTRT0V2 CASE F DTL= ',DTL,' DTH= ',
     &                 DTH,' DTFS1= ',DTFS1,' DXFS= ',DXFS1
               ENDIF
               IF (DXFS1.GT.BNDTOL) THEN
                  INTERSECT = .FALSE.
                  WRITE(6,*) 'PROBLEM CASE F DXFS1 = ',DXFS1,' WRONG'
                  CALL EXIT(1)
               ELSE
                  DTL = DTFS1
                  DTH = DT
                  DTFI= RTSAFE(RESDXF2,DTL,DTH,DXINCRTOL,DXNLTOL,
     &                 MAXIT,RPAR,IPAR)
                  DTF = DTFI
                  IF (IDEBUG.GT.0) THEN
C     MWF DEBUG
C UPDATE POSITION RELATIVE TO FACE
                     CALL DXUPDATEVXT(ZEROTOL,DTFI,V0F,VX0,VTF,VXT,DXFI)
                     DXFI = DXF0 + DXFI
                     WRITE(6,*)'STEPDTRT0V2 CASE F DTL= ',DTL,
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
               WRITE(6,*)'STEPDTRT0V2 INTERSECTED FACE= ',J,' DIR= ',
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
      
      DO I=1,NEQ
         CALL DXUPDATEVXT(ZEROTOL,DT,VTEMP(I),VX0,VT(I),VXT,DXOUT(I))
         XOUT(I) = XS(I) + DXOUT(I)
      ENDDO

      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'STEPDTRT0V2 AFTER FACE LOOP T= ',T,' DT= ',DT
         WRITE(6,*)'XS= ',(XS(I),I=1,NEQ)
         WRITE(6,*)'XOUT= ',(XOUT(I),I=1,NEQ)
      ENDIF

      RETURN
      END
C
C

C
C UPDATE FORMULAS ARE GIVEN BY EQUATIONS (12)--(15) IN RUSSELL, HEALY 00
C BUT ACTUAL IMPLEMENTATION GIVEN IN EQUATIONS (16)--(19)
      SUBROUTINE DXUPDATEVXT
     I     (ZEROTOL,DT,V0,VX,VT,VXT,
     O     DX)
C
      IMPLICIT NONE
C EPSILON, MACHINE_EPS
      DOUBLE PRECISION ZEROTOL
C \Delta t, v_0, v_x, v_t, v_{xt}
      DOUBLE PRECISION DT,V0,VX,VT,VXT
C x-x_0
      DOUBLE PRECISION DX
C LOCAL VARIABLES
C DELTA CUTOFF FOR v_xt AND S=v_x\Delta t FROM FORMULA 18 AND 19
      DOUBLE PRECISION DELTA,S,SMAG,SE
C     \frac{e^s-1}{v_x}, \frac{e^s-(1+s)}{v_x^2}
      DOUBLE PRECISION TERM1,TERM2
C     \frac{1}{v_x}\left(\frac{1}{2}\Delta t^2 e^s - \frac{e^s-(1+s)}{v^2_x}\right)
      DOUBLE PRECISION TERM3
C     \frac{e^s - (1+s+\frac{1}{2}s^2 + \frac{1}{6}s^3}{v_x^4}
      DOUBLE PRECISION TERM4
C     \frac{e^s - (1+s+\frac{1}{2}s^2)}{v^3_x}
      DOUBLE PRECISION TERM5
C     TERMS IN EVALUATING DELTA 16
C     |2\epsilon v_t \max[1,e^{v_x\Delta t}]|
      DOUBLE PRECISION DELTATERM1
C     \left|\frac{v_t v_x\Delta t^4}{8} \frac{e^s - (1+s+\frac{1}{2}s^2)}{v^3_x}\right|
      DOUBLE PRECISION DELTATERM2
C     \left|\frac{v_0\Delta t^4}{8}\frac{e^s-1}{v_x}\right | 
      DOUBLE PRECISION DELTATERM3
C     |\frac{v_t\Delta t^6}{48}|
C     
      DOUBLE PRECISION EPSLON,ONETHIRD,ONESIXTH,ONEEIGHTH,ONE24

C MACHINE EPSILON, SET FOR 64 BIT INTEL MAC RIGHT NOW
      EPSLON = 1.0D-16

      ONETHIRD = 1.D0/3.D0
      ONESIXTH = 1.D0/6.D0
      ONEEIGHTH= 1.D0/8.D0
      ONE24    = 1.D0/24.D0
      S    = VX*DT
      SMAG = DABS(S)
      SE   = DEXP(S)
C EVALUATE TERMS ON PAGE 6, RH00
      IF (SMAG.LT.(6.D0*ZEROTOL)**ONETHIRD) THEN
         TERM1 = DT*(1.D0 + 0.5*S)
      ELSE
         TERM1 = (SE-1.D0)/VX
      ENDIF
      IF (SMAG.LT.(24.D0*ZEROTOL)**(0.25D0)) THEN
         TERM2 = DT*DT*(0.5 + S*ONESIXTH)
      ELSE
         TERM2 = (SE-(1.D0+S))/(VX*VX)
      ENDIF
      IF (SMAG.LT.(40.D0*ONETHIRD*ZEROTOL)**(0.2D0)) THEN
         TERM3 = (DT**3)*(ONETHIRD + 5.D0*S*ONE24)
      ELSE
         TERM3 = (0.5*DT*DT*SE - (SE - (1.D0+S))/(VX*VX))/VX
      ENDIF
      IF (SMAG.LT.(720.D0*ZEROTOL)**ONESIXTH) THEN
         TERM4 = DT**(4.D0)*(ONE24 + S*ONESIXTH*0.05D0)
      ELSE
         TERM4 = (SE - (1.D0 + S + 0.5*S**2 + ONESIXTH*S**3))/(VX**4)
      ENDIF
      IF (SMAG.LT.(120.D0*ZEROTOL)**(0.2D0)) THEN
         TERM5 = DT**(3.D0)*(ONESIXTH + ONE24*S)
      ELSE
         TERM5 = (SE - (1.D0 + S + 0.5*S*S))/(VX**3)
      ENDIF

C EQUATION 19
      DELTATERM1 = DABS(2.D0*ZEROTOL*VT*MAX(1.D0,SE))
      DELTATERM2 = DABS(ONEEIGHTH*VT*VX*(DT**4)*TERM5)
      DELTATERM3 = DABS(ONEEIGHTH*V0*(DT**4)*TERM1) + 
     &     DABS(0.5*ONE24*VT*DT**(6))

      DELTA = DELTATERM1**(ONETHIRD)/((DELTATERM2+DELTATERM3)**ONETHIRD)

      IF (VXT.GT.DELTA) THEN
         CALL DXVTPOS(DT,V0,VX,VT,VXT,DX)
      ELSE IF (VXT.LT.-DELTA) THEN
         CALL DXVTNEG(DT,V0,VX,VT,VXT,DX)
C      ELSE IF (DABS(VX).LE.ZEROTOL) THEN
C         DX = V0*DT + 0.5*VT*DT*DT
      ELSE
         DX = V0*TERM1 + VT*TERM2 + VXT*(V0*TERM3 + 
     &        VT*(0.5*DT*DT*TERM2 - 3.D0*TERM4))

      ENDIF

      RETURN
      END

C EQUATION 12 FROM RUSSELL HEALY 2000
      SUBROUTINE DXVTPOS
     I     (DT,V0,VX,VT,VXT,
     O     DX)
      IMPLICIT NONE
C \Delta t, v_0, v_x, v_t, v_{xt}
      DOUBLE PRECISION DT,V0,VX,VT,VXT
C x-x_0
      DOUBLE PRECISION DX
C LOCAL VARIABLES
      DOUBLE PRECISION DPI,TMP1,TMP2,TMP3
C
      DPI = 3.14159265358979311599796

      IF (VXT.LE.0.D0) THEN
         WRITE(6,*) 'DXVTPOS REQUIRES VXT > 0, GOT ',VXT
         CALL EXIT(1)
      ENDIF
      TMP1 = VT/VXT*(DEXP(0.5D0*VXT*DT*DT + VX*DT)-1.D0)
      TMP2 = (V0-VX*VT/VXT)*DEXP(0.5D0*VXT*(DT + VX/VXT)**2)*
     &     DSQRT(0.5D0*DPI/VXT)
      TMP3 = DERFC(DSQRT(0.5D0*VXT)*VX/VXT) - 
     &     DERFC(DSQRT(0.5D0*VXT)*(DT + VX/VXT))
      
      DX = TMP1 + TMP2*TMP3

      RETURN
      END

C EQUATION 13 FROM RUSSELL HEALY 2000
      SUBROUTINE DXVTNEG
     I     (DT,V0,VX,VT,VXT,
     O     DX)
      IMPLICIT NONE
C \Delta t, v_0, v_x, v_t, v_{xt}
      DOUBLE PRECISION DT,V0,VX,VT,VXT
C x-x_0
      DOUBLE PRECISION DX
C LOCAL VARIABLES
      DOUBLE PRECISION TMP1,TMP2,TMP3,VXTMAG,
     &     LOW,HIG,LOWINT,HIGINT
C CALCULATES int_{0}^{x} e^{s^2} \ds 
      DOUBLE PRECISION DINTE2

      IF (VXT.GE.0.D0) THEN
         WRITE(6,*) 'DXVTPOS REQUIRES VXT < 0, GOT ',VXT
         CALL EXIT(1)
      ENDIF
      VXTMAG = DABS(VXT)
      TMP1   = VT/VXT*(DEXP(0.5D0*VXT*DT*DT + VX*DT)-1.D0)
      TMP2   = (V0-VX*VT/VXT)*DEXP(0.5D0*VXT*(DT + VX/VXT)**2)*
     &     DSQRT(2.D0/VXTMAG)
      LOW    = DSQRT(0.5D0*VXTMAG)*VX/VXT
      HIG    = DSQRT(0.5D0*VXTMAG)*(DT + VX/VXT)
      LOWINT = DINTE2(LOW)
      HIGINT = DINTE2(HIG)
      TMP3   = HIGINT-LOWINT
      
      DX     = TMP1 + TMP2*TMP3
      
      RETURN 
      END

C
C
C RESIDUAL ROUTINE FOR INTERSECTION WITH FACE f
C WITH BILINEAR (X-T) VARIATION
C  (X(dt) - X_f).N_f = R
C  dx_0f + dx_f(dt)  = R
C
C  dx_0f = (X_0-X_f).N_f
C  dx_f(dt)  = somewhat complicated function in eqns (12)--(15) RH00
C
C  Jacobian, R' is just d(dx_f)/d(dt) = v_f  
      SUBROUTINE RESDXF2
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
C XPAR(6) = v_xt
C LOCAL VARIABLES
      DOUBLE PRECISION DX0F,DXF,V0F,VTF,VX,VXT,VF
      DOUBLE PRECISION ZEROTOL
      
      DX0F= XPAR(1)
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)
      VXT     = XPAR(6)
      
      CALL DXUPDATEVXT(ZEROTOL,DT,V0F,VX,VTF,VXT,DXF)
      VF = V0F + VX*DXF + DT*VTF + VXT*DXF*DT

      R = DX0F + DXF
      DR= VF

      RETURN
      END
C
C RESIDUAL ROUTINE FOR VELOCITY CRITICAL POINT
C  V(X(dt),dt).N_f = R
C  [V_0 + v_x (X(dt)-X_0) + V_t dt + v_xt*(X(dt)-X_0)*dt].N_f      = R 
C  v_f(dt) = v_0f + v_x dx_f(dt) + v_tf dt  + v_xtf*dx_f(dt)*dt    = R 
C  dx_f(dt)= somewhat complicated function given by equations (12)--(15) in RH00
C  v_0f = V_0 .N_f, v_tf = V_t . N_f
C  
C DERIVATIVE IS
C v_f'(dt) = v_tf + v_x*dx_f'(dt) + v_xt*dx_f(dt) + v_xt*dx_f'*dt 
C          = v_tf + v_xt*dx_f + (v_x + v_xt*dt)*v_f(dt) 
C
C
C 
      SUBROUTINE RESVF2
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
C XPAR(6) = v_xt
C LOCAL VARIABLES
      DOUBLE PRECISION VF,DXF,V0F,VTF,VX,VXT
      DOUBLE PRECISION ZEROTOL
            
      V0F = XPAR(2)
      VTF = XPAR(3)
      VX  = XPAR(4)
      ZEROTOL = XPAR(5)
      VXT     = XPAR(6)

      CALL DXUPDATEVXT(ZEROTOL,DT,V0F,VX,VTF,VXT,DXF)

      VF  = V0F + VX*DXF + VTF*DT + VXT*DT*DXF
      R   = VF
      DR = VTF + VXT*DXF + (VX + VXT*DT)*VF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C AUXILIARY ROUTINES
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
C input
C 0 converge on |dx| or |f|
C 1 converge on both |f| and |dx| 
C 3 
C output
C -1 did not converge
C  0 converge on |dx| only
C  1 converge on |dx| and |facc| 
C  2 converge on |f| only
      ipar(3) = -1
C number of iterations
      ipar(2) = 0
      call funcd(x1,fl,df,xpar,ipar)
      call funcd(x2,fh,df,xpar,ipar)
       
      if((fl.gt.facc.and.fh.gt.facc).or.
     &     (fl.lt.(-facc).and.fh.lt.(-facc))) then
         write(6,*) 'root must  be bracketed in rtsafe'
         write(6,*) 'x1= ',x1,' x2= ',x2, 
     &        ' fl= ',fl,' fh= ',fh,' facc= ',facc 
         call exit(1)
      endif

      if(dabs(fl).le.facc) then
         rtsafe = x1
         ipar(3) = 2
         return
      else if(dabs(fh).le.facc) then
         rtsafe = x2
         ipar(3)= 2
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
               ipar(3) = 0
C shouldn't be true
               if (dabs(fl).le.facc) ipar(3)=1
               return
            endif
         else
C Newton step acceptable
            dxold=dx
            dx=f/df
            temp=rtsafe
            rtsafe=rtsafe-dx
            if(temp.eq.rtsafe) then
               ipar(3) = 0
               if (dabs(f).le.facc) ipar(3)=1
               return
            endif
         endif
C 
         call funcd(rtsafe,f,df,xpar,ipar)
C convergence test, move after feval to use both tests
         if (dabs(dx).lt.xacc.and.dabs(f).le.facc) then
            ipar(3) = 1
            return
         else if (ipar(1).eq.0.and.dabs(dx).lt.xacc) then
            ipar(3) = 0
            return
         else if (ipar(1).eq.0.and.dabs(f).le.facc) then
            ipar(3) = 2
            return
         endif
C maintain root bracketing
         if(f.lt.0.D0) then
            xl=rtsafe
         else
            xh=rtsafe
         endif
      enddo
C      if (dabs(dx).lt.xacc.or.dabs(f).le.facc) then
C         write(6,*)'rtsafe exceeded maximum its= ',maxit, 
C     &        'dx= ',dx,' f= ',f, 'not quitting!'
C         ipar(3) = -1
C      else
      write(6,*)'rtsafe exceeded maximum its= ',maxit, 
     &     'dx= ',dx,' f= ',f, ' quitting!'
      call exit(1)
C      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C UTILITY ROUTINE FOR EVALUATING UPDATE FORMULAS IN BILINEAR (X-T) CASE
      
      FUNCTION DAW(XX)
C----------------------------------------------------------------------
C
C This function program evaluates Dawson's integral, 
C
C                       2  / x   2
C                     -x   |    t
C             F(x) = e     |   e    dt
C                          |
C                          / 0
C
C   for a real argument x.
C
C   The calling sequence for this function is 
C
C                   Y=DAW(X)
C
C   The main computation uses rational Chebyshev approximations
C   published in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek
C   and Thacher.  This transportable program is patterned after the
C   machine-dependent FUNPACK program DDAW(X), but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   approximations that are theoretically accurate to about 19
C   significant decimal digits.  The accuracy achieved depends on the
C   arithmetic system, the compiler, the intrinsic functions, and
C   proper selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XINF   = largest positive machine number
C   XMIN   = the smallest positive machine number.
C   EPS    = smallest positive number such that 1+eps > 1.
C            Approximately  beta**(-p), where beta is the machine
C            radix and p is the number of significant base-beta
C            digits in a floating-point number.
C   XMAX   = absolute argument beyond which DAW(X) underflows.
C            XMAX = min(0.5/xmin, xinf).
C   XSMALL = absolute argument below DAW(X)  may be represented
C            by X.  We recommend XSMALL = sqrt(eps).
C   XLARGE = argument beyond which DAW(X) may be represented by
C            1/(2x).  We recommend XLARGE = 1/sqrt(eps).
C
C     Approximate values for some important machines are
C
C                        beta  p     eps     xmin       xinf  
C
C  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
C  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
C  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
C  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
C                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
C   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307
C
C                         XSMALL     XLARGE     XMAX    
C
C  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
C  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
C  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
C  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
C                (D.P.)  3.73E-09   2.68E+08  1.70E+38
C   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns 0.0 for |X| > XMAX.
C
C Intrinsic functions required are:
C
C     ABS
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: June 15, 1988
C
C----------------------------------------------------------------------
      INTEGER I
CS    REAL
      DOUBLE PRECISION
     1     DAW,FRAC,HALF,ONE,ONE225,P1,P2,P3,P4,Q1,Q2,Q3,Q4,SIX25,
     2     SUMP,SUMQ,TWO5,W2,X,XX,Y,XLARGE,XMAX,XSMALL,ZERO
      DIMENSION P1(10),P2(10),P3(10),P4(10),Q1(10),Q2(9),Q3(9),Q4(9)
C----------------------------------------------------------------------
C  Mathematical constants.
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE/0.0E0,0.5E0,1.0E0/,
CS   1     SIX25,ONE225,TWO5/6.25E0,12.25E0,25.0E0/
      DATA ZERO,HALF,ONE/0.0D0,0.5D0,1.0D0/,
     1     SIX25,ONE225,TWO5/6.25D0,12.25D0,25.0D0/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XSMALL/2.44E-04/, XLARGE/4.10E+03/, XMAX/4.25E+37/
      DATA XSMALL/1.05D-08/, XLARGE/9.49D+07/, XMAX/2.24D+307/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation for  |x| < 2.5
C----------------------------------------------------------------------
CS    DATA P1/-2.69020398788704782410E-12, 4.18572065374337710778E-10,
CS   1        -1.34848304455939419963E-08, 9.28264872583444852976E-07,
CS   2        -1.23877783329049120592E-05, 4.07205792429155826266E-04,
CS   3        -2.84388121441008500446E-03, 4.70139022887204722217E-02,
CS   4        -1.38868086253931995101E-01, 1.00000000000000000004E+00/
CS    DATA Q1/ 1.71257170854690554214E-10, 1.19266846372297253797E-08,
CS   1         4.32287827678631772231E-07, 1.03867633767414421898E-05,
CS   2         1.78910965284246249340E-04, 2.26061077235076703171E-03,
CS   3         2.07422774641447644725E-02, 1.32212955897210128811E-01,
CS   4         5.27798580412734677256E-01, 1.00000000000000000000E+00/
      DATA P1/-2.69020398788704782410D-12, 4.18572065374337710778D-10,
     1     -1.34848304455939419963D-08, 9.28264872583444852976D-07,
     2        -1.23877783329049120592D-05, 4.07205792429155826266D-04,
     3        -2.84388121441008500446D-03, 4.70139022887204722217D-02,
     4        -1.38868086253931995101D-01, 1.00000000000000000004D+00/
      DATA Q1/ 1.71257170854690554214D-10, 1.19266846372297253797D-08,
     1         4.32287827678631772231D-07, 1.03867633767414421898D-05,
     2         1.78910965284246249340D-04, 2.26061077235076703171D-03,
     3         2.07422774641447644725D-02, 1.32212955897210128811D-01,
     4         5.27798580412734677256D-01, 1.00000000000000000000D+00/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  x in [2.5, 3.5)
C----------------------------------------------------------------------
CS    DATA P2/-1.70953804700855494930E+00,-3.79258977271042880786E+01,
CS   1         2.61935631268825992835E+01, 1.25808703738951251885E+01,
CS   2        -2.27571829525075891337E+01, 4.56604250725163310122E+00,
CS   3        -7.33080089896402870750E+00, 4.65842087940015295573E+01,
CS   4        -1.73717177843672791149E+01, 5.00260183622027967838E-01/
CS    DATA Q2/ 1.82180093313514478378E+00, 1.10067081034515532891E+03,
CS   1        -7.08465686676573000364E+00, 4.53642111102577727153E+02,
CS   2         4.06209742218935689922E+01, 3.02890110610122663923E+02,
CS   3         1.70641269745236227356E+02, 9.51190923960381458747E+02,
CS   4         2.06522691539642105009E-01/
      DATA P2/-1.70953804700855494930D+00,-3.79258977271042880786D+01,
     1         2.61935631268825992835D+01, 1.25808703738951251885D+01,
     2        -2.27571829525075891337D+01, 4.56604250725163310122D+00,
     3        -7.33080089896402870750D+00, 4.65842087940015295573D+01,
     4        -1.73717177843672791149D+01, 5.00260183622027967838D-01/
      DATA Q2/ 1.82180093313514478378D+00, 1.10067081034515532891D+03,
     1        -7.08465686676573000364D+00, 4.53642111102577727153D+02,
     2         4.06209742218935689922D+01, 3.02890110610122663923D+02,
     3         1.70641269745236227356D+02, 9.51190923960381458747D+02,
     4         2.06522691539642105009D-01/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  x in [3.5, 5.0]
C----------------------------------------------------------------------
CS    DATA P3/-4.55169503255094815112E+00,-1.86647123338493852582E+01,
CS   1        -7.36315669126830526754E+00,-6.68407240337696756838E+01,
CS   2         4.84507265081491452130E+01, 2.69790586735467649969E+01,
CS   3        -3.35044149820592449072E+01, 7.50964459838919612289E+00,
CS   4        -1.48432341823343965307E+00, 4.99999810924858824981E-01/
CS    DATA Q3/ 4.47820908025971749852E+01, 9.98607198039452081913E+01,
CS   1         1.40238373126149385228E+01, 3.48817758822286353588E+03,
CS   2        -9.18871385293215873406E+00, 1.24018500009917163023E+03,
CS   3        -6.88024952504512254535E+01,-2.31251575385145143070E+00,
CS   4         2.50041492369922381761E-01/
      DATA P3/-4.55169503255094815112D+00,-1.86647123338493852582D+01,
     1        -7.36315669126830526754D+00,-6.68407240337696756838D+01,
     2         4.84507265081491452130D+01, 2.69790586735467649969D+01,
     3        -3.35044149820592449072D+01, 7.50964459838919612289D+00,
     4        -1.48432341823343965307D+00, 4.99999810924858824981D-01/
      DATA Q3/ 4.47820908025971749852D+01, 9.98607198039452081913D+01,
     1         1.40238373126149385228D+01, 3.48817758822286353588D+03,
     2        -9.18871385293215873406D+00, 1.24018500009917163023D+03,
     3        -6.88024952504512254535D+01,-2.31251575385145143070D+00,
     4         2.50041492369922381761D-01/
C----------------------------------------------------------------------
C  Coefficients for R(9,9) approximation in J-fraction form
C     for  |x| > 5.0
C----------------------------------------------------------------------
CS    DATA P4/-8.11753647558432685797E+00,-3.84043882477454453430E+01,
CS   1        -2.23787669028751886675E+01,-2.88301992467056105854E+01,
CS   2        -5.99085540418222002197E+00,-1.13867365736066102577E+01,
CS   3        -6.52828727526980741590E+00,-4.50002293000355585708E+00,
CS   4        -2.50000000088955834952E+00, 5.00000000000000488400E-01/
CS    DATA Q4/ 2.69382300417238816428E+02, 5.04198958742465752861E+01,
CS   1         6.11539671480115846173E+01, 2.08210246935564547889E+02,
CS   2         1.97325365692316183531E+01,-1.22097010558934838708E+01,
CS   3        -6.99732735041547247161E+00,-2.49999970104184464568E+00,
CS   4         7.49999999999027092188E-01/
      DATA P4/-8.11753647558432685797D+00,-3.84043882477454453430D+01,
     1        -2.23787669028751886675D+01,-2.88301992467056105854D+01,
     2        -5.99085540418222002197D+00,-1.13867365736066102577D+01,
     3        -6.52828727526980741590D+00,-4.50002293000355585708D+00,
     4        -2.50000000088955834952D+00, 5.00000000000000488400D-01/
      DATA Q4/ 2.69382300417238816428D+02, 5.04198958742465752861D+01,
     1         6.11539671480115846173D+01, 2.08210246935564547889D+02,
     2         1.97325365692316183531D+01,-1.22097010558934838708D+01,
     3        -6.99732735041547247161D+00,-2.49999970104184464568D+00,
     4         7.49999999999027092188D-01/
C----------------------------------------------------------------------
      X = XX
      IF (ABS(X) .GT. XLARGE) THEN
            IF (ABS(X) .LE. XMAX) THEN
                  DAW = HALF / X
               ELSE
                  DAW = ZERO
            END IF
         ELSE IF (ABS(X) .LT. XSMALL) THEN
            DAW = X
         ELSE
            Y = X * X
            IF (Y .LT. SIX25) THEN
C----------------------------------------------------------------------
C  ABS(X) .LT. 2.5 
C----------------------------------------------------------------------
                  SUMP = P1(1)
                  SUMQ = Q1(1)
                  DO 100 I = 2, 10
                     SUMP = SUMP * Y + P1(I)
                     SUMQ = SUMQ * Y + Q1(I)
  100             CONTINUE
                  DAW = X * SUMP / SUMQ
               ELSE IF (Y .LT. ONE225) THEN
C----------------------------------------------------------------------
C  2.5 .LE. ABS(X) .LT. 3.5 
C----------------------------------------------------------------------
                  FRAC = ZERO
                  DO 200 I = 1, 9
  200                FRAC = Q2(I) / (P2(I) + Y + FRAC)
                  DAW = (P2(10) + FRAC) / X
               ELSE IF (Y .LT. TWO5) THEN
C----------------------------------------------------------------------
C  3.5 .LE. ABS(X) .LT. 5.0 
C---------------------------------------------------------------------
                  FRAC = ZERO
                  DO 300 I = 1, 9
  300                FRAC = Q3(I) / (P3(I) + Y + FRAC)
                  DAW = (P3(10) + FRAC) / X
               ELSE
C----------------------------------------------------------------------
C  5.0 .LE. ABS(X) .LE. XLARGE 
C------------------------------------------------------------------
                  W2 = ONE / X / X
                  FRAC = ZERO
                  DO 400 I = 1, 9
  400                FRAC = Q4(I) / (P4(I) + Y + FRAC)
                  FRAC = P4(10) + FRAC
                  DAW = (HALF + HALF * W2 * FRAC) / X
            END IF
      END IF
      RETURN
C---------- Last line of DAW ----------
      END


      FUNCTION DINTE2(XX)
      IMPLICIT NONE
C WRAPPER FOR FUNCTION
C 
C  F(X) = \int_{0}^{x}e^{t^2}dt
C
C USING DAWSON'S INTEGRAL
C  
C  D(X) = e^{-x^2}\int_{0}^{x}e^{t^2}dt
C
      DOUBLE PRECISION DAW
      DOUBLE PRECISION DINTE2,XX,DAWVAL,DTMP
      
      DTMP = DEXP(XX*XX)
      DAWVAL = DAW(XX)
      DINTE2 = DAWVAL*DTMP

      RETURN
      END

C ======================================================================
C CONVERT TO PEARCE'S LATEST VERSION
C ======================================================================
C
C 
C
      SUBROUTINE ELTRAK123ANEW
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
C
      DIMENSION XI(3),DI(3),XI_S(3),DI_S(3)
      DIMENSION IXI(3),IDI(4),IXI_S(3),IDI_S(4)
      
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
CMWF HACK TEMPORARILY COMPUTE DL HERE
      CALL INTRP123A
     I     (MAXEQ,MAXND,NEQ,NODE, XS, XW,
     O     DN_S,DJAC)
      DL=DABS(DJAC)**(DEQ)

CMWF NEED TO COMPUTE DN_S IDI_S,IXI_S,XI_S,DI_S FOR PHI CALCS
      CALL INTRP123
     I      (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ZEROTOL,DL, XS, XW,
     O       DN_S,IADJUST,XI_S,DI_S,IXI_S,IDI_S)
CMWF IGNORE IADJST
      IF (IDEBUG.GT.0) THEN
CMWF DEBUG
         WRITE(6,*)'ELTRAK123ANEW CALLING STEP IVFLAG= ',IVFLAG,
     &        ' DT= ',DT,' DL= ',DL,' XS= ',(XS(I),I=1,NEQ)
         WRITE(6,*)' DN_S= ',(DN_S(I),I=1,NODE)
         DO J=1,NODE
            WRITE(6,*)' ',J,' XW= ',(XW(II,J),II=1,NEQ)
         ENDDO
      ENDIF

      IF (IVFLAG.EQ.1) THEN
C RT0, LINEAR DEPENDENCE IN TIME
         CALL STEPDTRT0V1
     I        (MAXEQ,MAXND,NEQ,NODE,
     I        DEQ, DN_SAFE,
     I        ZEROTOL,DIR,T1,T2,
     I        XW,VT1W,VT2W,
     I        FNORMALS,XFBAR,
     M        T,DT,XS,XOUT5,XTEMP,
     M        VTEMP,DN_S,DN,DDN)
      ELSE IF (IVFLAG.EQ.2) THEN
C RT0, BILINEAR IN X-T
         CALL STEPDTRT0V2
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

      PHI=0.0D0
      CALL INTRP123
     I      (MAXEQ,MAXND,NEQ,NODE,DN_SAFE,ZEROTOL,DL, XOUT5, XW,
     O       DN,IADJUST,XI,DI,IXI,IDI)
C
C < NOTE > ADJUST XOUT5 WHEN NECESSARY (I.E., IADJUST = 1)
C
      IF (IDEBUG.GT.2) THEN 
         WRITE(6,*)'AFTER STEP INTRP123 IADJUST= ',IADJUST,
     &        ' DL= ',DL
         WRITE(6,*)'DN_S= ',(DN_S(I),I=1,NODE)
         WRITE(6,*)' XOUT5= ',(XOUT5(I),I=1,NEQ)
         WRITE(6,*)'DN= ',(DN(I),I=1,NODE)
      ENDIF
      IF(IADJUST.EQ.1)THEN
         DO I=1,NEQ
            XOUT5(I)=0.0E0
            DO J=1,NODE
               XOUT5(I)=XOUT5(I)+DN(J)*XW(I,J)
            ENDDO
         ENDDO
      ENDIF                 

C
C === COMPUTE PHI
C
      CALL PHI_COMP
     I      (MAXND,NODE,NEQ, 
     I       DN_SAFE, 
     I       DN_S,DN,XI,XI_S,DI,DI_S, IXI_S,IDI_S,
     O       IDSDT,I1,I2,I3,PHI)
CMWF
      IF(IDEBUG.GT.2) THEN
         WRITE(6,*)'AFTER STEP PHI_COMP IDSDT= ',IDSDT,
     &        ' PHI= ',PHI,' I1= ',I1,' I2= ',I2,' I3 ',I3
      ENDIF

CMWF NO TRACKING THROUGH ELEMENT
      IF(IDSDT.EQ.0)RETURN
C
C === IF PHI IS GREATER THAN 1 ==> REDUCE TIMESTEP
C
      IF(PHI.GT.1.0E0)THEN
         WRITE(6,*)' PROBLEM IN ELTRAK123ANEW IDSDT= ',IDSDT,
     &        ' PHI= ',PHI, 'SHOULDNT HAPPEN'
         WRITE(6,*)' DT= ',DT,' XS= ',(XS(I),I=1,NEQ)
         WRITE(6,*)' XOUT5= ',(XOUT5(I),I=1,NEQ)
         WRITE(6,*)' DN= ',(DN(I),I=1,NODE)
         CALL EXIT(1)
      ENDIF
C
C B. WHEN THE ENDING LOCATION IS EITHER WITHIN THE ELEMENT OR 
C    ON THE ELEMENT BOUDNARY
C ==> UPDATE INFORMATION FOR THE SUCCESSIVE PT
C
      T=T+DT
      TPT(IPT)=T
      DO I=1,NEQ
         XS(I)=XOUT5(I)
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
      CALL EB_CHECK
     I     (NEQ,NODE, IDI,IXI, 
     O     I1,I2,I3)
      IF(I1.NE.0) THEN 
         IF (IDEBUG.GT.2) THEN
            WRITE(6,*)' ELTRAK123ANEW LEAVING ON BOUNDARY IDSDT= ',
     &             IDSDT,' I1= ',I1,' I2= ',I2,' I3= ',I3
         ENDIF
         RETURN
      ENDIF

C MWF DEBUG
      IF(IDSDT.EQ.-1) THEN
         IF (IDEBUG.GT.0) THEN
C MWF DEBUG
            WRITE(6,*)' ELTRAK123A INSIDE ELE RETURN IDSDT= ',IDSDT, 
     &           ' I1= ',I1,' I2= ',I2,' I3= ',I3, 
     &           ' DN= ',(DN(II),II=1,NODE)
         ENDIF
         RETURN
      ENDIF
C MWF END DEBUG
C MWF THIS SHOULD NOT HAPPEN WITH ANALYTICAL TRACKING?
      WRITE(6,*)'PROBLEM ELTRAK123ANEW END WITHOUT CONCLUSION= ',
     &     ' T= ',T,' DT= ',DT
      WRITE(6,*)'XS= ',(XS(II),II=1,NEQ)
      WRITE(6,*)'XOUT5= ',(XOUT5(II),II=1,NEQ)
      CALL EXIT(1)
      
C 
C  999 CONTINUE
      RETURN
      END

