"""
This file declares a list of C++ objects from Chrono that can be used in
 cython files.
"""


from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)

cdef extern from "ProtChMoorings.h":

    # ------- CORE ------- #

    cdef cppclass ChVector3d:
        double x()
        double y()
        double z()
        ChVector3d(double x,
                 double y,
                 double z)
        ChVector3d()

    cdef cppclass ChVector3i:
        int x()
        int y()
        int z()
        ChVector3d(int x,
                 int y,
                 int z)
        ChVector3i()

    cdef cppclass ChWrenchd:
      ChVector3d force
      ChVector3d torque
      ChWrenchd()

    cdef cppclass ChQuaternion[double]:
        double e0()
        double e1()
        double e2()
        double e3()
        ChQuaternion(double e0,
                     double e1,
                     double e2,
                     double e3)
        ChQuaternion()

    cdef cppclass ChCoordsys[double]:
        ChVector3d pos
        ChQuaternion rot
        ChCoordsys()
        ChCoordsys(ChVector3d &mv,
                   ChQuaternion &mq)
        ChCoordsys(ChVector3d &mv,
                   double Alpha,
                   ChVector3d &mu)

    cdef cppclass ChMatrix33[double]:
        ChVector3d GetAxisX()
        ChVector3d GetAxisY()
        ChVector3d GetAxisZ()
        # void CopyFromMatrixT(ChMatrix matra)
        # ChMatrix33()
        # ChMatrix33(ChQuaternion quat)

    cdef cppclass ChMatrixDynamic[double]:
        ChMatrixDynamic()
        ChMatrixDynamic(const int row,
                        const int col)
        # ChMatrixDynamic operator=(const ChMatrix& matbis)
        # ChMatrixDynamic operator+(const ChMatrix& matbis)
        ChVector3d GetAxisX()
        ChVector3d GetAxisY()
        ChVector3d GetAxisZ()

    cdef cppclass ChTriangle:
        ChTriangle(const ChVector3d &mp1,
                   const ChVector3d &mp2,
                   const ChVector3d &mp3)

    cdef cppclass ChTriangleMesh:
        void addTriangle(const ChVector3d &vertex0,
                         const ChVector3d &vertex1,
                         const ChVector3d &vertex2)
        void addTriangle(const ChTriangle &atriangle)

    cdef cppclass ChTriangleMeshConnected(ChTriangleMesh):
        ChTriangleMeshConnected()
        vector[ChVector3d]& GetCoordsVertices()
        vector[ChVector3d]& GetCoordsNormals()
        vector[ChVector3i]& GetIndicesVertexes()
        void LoadWavefrontMesh(string filename,
                               bool load_normals=True,
                               bool load_uv=False)

    cdef cppclass ChCollisionModel:
        bool AddTriangleMesh(shared_ptr[ChTriangleMesh] trimesh,
                             bool is_static,
                             bool is_convex,
                             const ChVector3d &pos,
                             const ChMatrix33 &rot,
                             double sphereswept_thickness)
        void SetEnvelope(double amargin)
        void SetSafeMargin(double amargin)
        void ClearModel()
        void BuildModel()

    ChQuaternion Q_from_AngAxis(double angle,
                                const ChVector3d &axis)

    # ------- PHYSICS ------- #

    cdef cppclass ChSystem:
        void Add(shared_ptr[ChPhysicsItem] newitem)
        void AddBody(shared_ptr[ChBody] newbody)
        double GetChTime()
        void SetupInitial()

    cdef cppclass ChSystemSMC:
        void Add(shared_ptr[ChPhysicsItem] newitem)
        void AddBody(shared_ptr[ChBody] newbody)
        double GetStep()

    cdef cppclass ChPhysicsItem:
        ChPhysicsItem()

    cdef cppclass ChFrame[double]:
        void SetPos(const ChVector3d& mpos) except +
        ChVector3d& GetPos()
        ChQuaternion& GetRot()
        void SetRot(ChQuaternion &rot) except +
        void SetPos(ChVector3d &pos) except +
        ChMatrix33& GetRotMat()
        ChVector3d GetRotAxis()
        double GetRotAngle()

    cdef cppclass ChFrameMoving[double](ChFrame):
        ChVector3d& GetPosDt()
        void SetPosDt(ChVector3d &pos_dt)
        ChQuaternion& GetRotDt()
        void SetRotDt(ChQuaternion &mrot_dt)
        ChVector3d& GetPosDt2()
        void SetPosDt2(ChVector3d &mpos_dtdt)
        ChQuaternion& GetRotDt2()
        void SetRotDt2(ChQuaternion &mrot_dtdt)
        ChMatrix33 GetRotMatDt()
        ChMatrix33 GetRotMatDt2()
        ChVector3d GetAngVelLocal()
        ChVector3d GetAngAccLocal()

    cdef cppclass ChBodyFrame(ChFrameMoving):
        ChBodyFrame()

    cdef cppclass ChBody(ChPhysicsItem, ChBodyFrame):
        ChBody() except +
        # void SetRot(ChQuaternion &rot) except +
        # void SetInertiaXX(ChVector3d &iner)
        # void SetInertiaXY(ChVector3d &iner)
        # const ChMatrix33& GetInertia()
        void SetBodyFixed(bool state) except +
        void SetMaterialSurface(const shared_ptr[ChMaterialSurface] &mnewsurf) except +
        shared_ptr[ChCollisionModel] GetCollisionModel()
        void EnableCollision(bool state)
        bool IsCollisionEnabled()
        # void SetMass(double newmass)
        # double GetMass()

    cdef cppclass ChBodyEasyBox(ChBody):
        ChBodyEasyBox(double Xsize,
                      double Ysize,
                      double Zsize,
                      double mdensity,
                      bool collide=False,
                      bool visual_asset=True)


    # ------- NODES ------- #

    cdef cppclass ChNodeBase:
        ChNodeBase()

    cdef cppclass ChNodeXYZ(ChNodeBase):
        void SetPos (const ChVector3d &mpos)
        const ChVector3d& GetPos()
        void SetPosDt (const ChVector3d &mposdt)
        const ChVector3d& GetPosDt()
        void SetPosDt2 (const ChVector3d &mposdtdt)
        const ChVector3d& GetPosDt2()

    cdef cppclass ChNodeFEAbase(ChNodeBase):
        ChNodeFEAbase()
        void Relax()
        void SetNoSpeedNoAcceleration()
        void SetFixed(bool mev)
        bool GetFixed()
        void SetIndex()
        int GetIndex()

    cdef cppclass ChNodeFEAxyz(ChNodeXYZ, ChNodeFEAbase):
        ChNodeFEAxyz()
        double GetMass()
        void SetMass(double mm)
        void SetForce(ChVector3d mf)
        ChVector3d& GetForce()

    cdef cppclass ChNodeFEAxyzD(ChNodeFEAxyz):
        const ChVector3d& GetSlope1()

    cdef cppclass ChNodeFEAxyzDD(ChNodeFEAxyzD):
        const ChVector3d& GetSlope2()

    cdef cppclass ChNodeFEAxyzrot(ChNodeFEAbase, ChBodyFrame):
        void SetForce(ChVector3d mf)
        ChVector3d& GetForce()
        void SetTorque(ChVector3d mf)
        ChVector3d& GetTorque()

    cdef cppclass ChNodeFEAxyzP(ChNodeFEAbase)

    cdef cppclass ChNodeFEAcurv(ChNodeFEAbase)


    # ------- ELEMENTS ------- #
    cdef cppclass ChElementBase:
        shared_ptr[ChNodeFEAbase] GetNodeN(int n)

    cdef cppclass ChElementGeneric(ChElementBase):
        ChElementGeneric()

    cdef cppclass ChElementBeam(ChElementGeneric):
        # void EvaluateSectionDisplacement(const double eta,
        #                                  const ChMatrix &displ,
        #                                  ChVector3d &u_displ,
        #                                  ChVector3d &u_rotaz)
        # void EvaluateSectionFrame(const double eta,
        #                           const ChMatrix &displ,
        #                           ChVector3d &u_displ,
        #                           ChQuaternion &rot)
        # void EvaluateSectionForceTorque(const double eta,
        #                                 const ChMatrix &displ,
        #                                 ChVector3d &Fforce,
        #                                 ChVector3d &Mtorque)
        # void EvaluateSectionStrain(const double eta,
        #                            const ChMatrix &displ,
        #                            ChVector3d &StrainV)
        double GetMass()
        double GetRestLength()
        void SetRestLength(double ml)

    cdef cppclass ChElementCableANCF(ChElementBeam):
        void SetNodes(shared_ptr[ChNodeFEAxyzD] nodeA, shared_ptr[ChNodeFEAxyzD] nodeB)
        # void SetSection()
        void SetNodes(shared_ptr[ChNodeFEAxyzD] nodeA, shared_ptr[ChNodeFEAxyzD] nodeB)

    cdef cppclass ChElementBeamANCF:
        void SetNodes(shared_ptr[ChNodeFEAxyzDD] nodeA, shared_ptr[ChNodeFEAxyzDD] nodeB, shared_ptr[ChNodeFEAxyzDD] nodeC)
        void SetDimensions(double lenX,
                           double beam_h,
                           double beam_w)
        shared_ptr[ChNodeFEAxyzDD] GetNodeA()
        shared_ptr[ChNodeFEAxyzDD] GetNodeB()
        shared_ptr[ChNodeFEAxyzDD] GetNodeC()
        void setAlphaDamp(double a)
        double GetLengthX() const
        double GetMass()

    cdef cppclass ChElementBeamEuler:
        ChElementBeamEuler()

    # ------- SECTIONS ------- #
    cdef cppclass ChBeamSection:
        bool IsCircular()
        void SetCircular(bool ic)


    cdef cppclass ChBeamSectionCable(ChBeamSection):
        ChBeamSectionCable() except +
        void SetArea(const double ma)
        double GetArea() const
        void SetI(double ma)
        double GetI() const
        void SetDiameter(double diameter)
        void SetDensity(double md)
        double GetDensity() const
        void SetYoungModulus(double mE)
        double GetYoungModulus() const
        void SetBeamRaleyghDamping(double mr)
        double GetBeamRaleyghDamping()



    cdef cppclass ChMesh:
        void SetAutomaticGravity(bool mg,
                                 int num_points=1)

    cdef cppclass ChMaterialSurface:
        ChMaterialSurface() except +

    cdef cppclass ChContactMaterial:
        ChContactMaterial() except +

    cdef cppclass ChContactMaterialSMC(ChContactMaterial):
        void SetYoungModulus(float val)
        void SetPoissonRatio(float val)
        void SetSfriction(float val)
        void SetKfriction(float val)
        void SetFriction(float val)
        void SetRestitution(float val)
        void SetAdhesion(float val)
        void SetAdhesionMultDMT(float val)
        void SetKn(float val)
        void SetKt(float val)
        void SetGn(float val)
        void SetGt(float val)
        double GetYoungModulus()
        double GetPoissonRatio()
        double GetSfriction()
        double GetKfriction()
        double GetRestitution()
        double GetAdhesion()
        double GetAdhesionMultDMT()
        double GetKn()
        double GetKt()
        double GetGn()
        double GetGt()

    cdef cppclass ChContactSurface:
        ChMesh* GetMesh()

    cdef cppclass ChContactSurfaceNodeCloud(ChContactSurface):
        ChContactSurfaceNodeCloud()
        void AddNode(shared_ptr[ChNodeFEAxyz] mnode,
                     const double point_radius=0.001)
        void AddNode(shared_ptr[ChNodeFEAxyzrot] mnode,
                     const double point_radius=0.001)
        void AddAllNodes(const double point_radius)

    cdef cppclass ChLinkBase:
        ChWrenchd GetReaction1()
        ChWrenchd GetReaction2()

    cdef cppclass ChLink(ChLinkBase)

    cdef cppclass ChLinkMate(ChLink)

    cdef cppclass ChLinkMateGeneric(ChLinkMate)

    cdef cppclass ChLinkPointFrame(ChLinkBase):
        int Initialize(shared_ptr[ChNodeFEAxyz] node, shared_ptr[ChBodyFrame] body, ChVector3d* pos)
        ChVector3d GetReactionOnNode()
        ChVector3d GetReactionOnBody()
        #virtual
        #int Initialize(shared_ptr[ChNodeFEAxyz] node, shared_ptr[ChBodyFrame] body, ChVector3d *pos=0)

    cdef cppclass ChLinkPointPoint(ChLinkBase):
        #virtual
        ChLinkPointPoint()
        int Initialize(shared_ptr[ChNodeFEAxyz] anodeA, shared_ptr[ChNodeFEAxyz] anodeB)


cdef extern from "ChBodyAddedMass.h":
    cdef cppclass ChBodyAddedMass(ChBody):
        ChBodyAddedMass() except +
        void SetMass(double newmass)
        void SetInertia(ChMatrix33& newXInertia)
        void SetInertiaXX(ChVector3d& newXInertia)
        void SetInertiaXY(ChVector3d& newXInertia)
        double GetMass()
        const ChMatrix33& GetInertia()
        ChVector3d GetInertiaXX()
        ChVector3d GetInertiaXY()
        void SetMfullmass(ChMatrixDynamic Mfullmass_in)
        void SetInvMfullmass(ChMatrixDynamic inv_Mfullmass_in)
    ChBodyAddedMass * newChBodyAddedMass()
