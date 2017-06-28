
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)

cdef extern from "ChMoorings.h":

    # ------- CORE ------- #

    cdef cppclass ChVector[double]:
        double x()
        double y()
        double z()
        ChVector(double x, double y, double z)
        ChVector()

    cdef cppclass ChQuaternion[double]:
        double e0()
        double e1()
        double e2()
        double e3()
        ChQuaternion(double e0, double e1, double e2, double e3)
        ChQuaternion()

    cdef cppclass ChMatrix

    cdef cppclass ChMatrix33[double]:
        ChVector Get_A_Xaxis()
        ChVector Get_A_Yaxis()
        ChVector Get_A_Zaxis()


    # ------- PHYSICS ------- #

    cdef cppclass ChSystem:
        void Add(shared_ptr[ChPhysicsItem] newitem)
        void AddBody(shared_ptr[ChBody] newbody)

    cdef cppclass ChSystemSMC:
        void Add(shared_ptr[ChPhysicsItem] newitem)
        void AddBody(shared_ptr[ChBody] newbody)
        double GetStep()
        double GetChTime()
        void SetupInitial()

    cdef cppclass ChPhysicsItem:
        ChPhysicsItem()

    cdef cppclass ChFrame[double]:
        void SetPos(const ChVector[double]& mpos) except +
        ChVector& GetPos()
        ChQuaternion& GetRot()
        void SetRot(ChQuaternion &rot) except +
        void SetPos(ChVector &pos) except +
        ChMatrix33& GetA()
        ChVector GetRotAxis()
        double GetRotAngle()

    cdef cppclass ChFrameMoving[double](ChFrame):
        ChVector& GetPos_dt()
        ChQuaternion& GetRot_dt()
        ChVector& GetPos_dtdt()
        ChQuaternion& GetRot_dtdt()
        ChMatrix33 GetA_dt()
        ChMatrix33 GetA_dtdt()
        ChVector GetWvel_loc()
        ChVector GetWacc_loc()

    cdef cppclass ChBodyFrame(ChFrameMoving):
        ChBodyFrame()

    cdef cppclass ChBody(ChPhysicsItem, ChBodyFrame):
        ChBody() except +
        # void SetRot(ChQuaternion &rot) except +
        void SetInertiaXX(ChVector &iner)
        void SetInertiaXY(ChVector &iner)
        void SetBodyFixed(bool state) except +
        void SetMaterialSurface(const shared_ptr[ChMaterialSurface] &mnewsurf) except +
        void SetMass(double newmass)
        double GetMass()

    cdef cppclass ChBodyEasyBox(ChBody):
        ChBodyEasyBox(double Xsize, double Ysize, double Zsize, double mdensity, bool collide=False, bool visual_asset=True)


    # ------- NODES ------- #

    cdef cppclass ChNodeBase:
        ChNodeBase()

    cdef cppclass ChNodeXYZ(ChNodeBase):
        void SetPos (const ChVector &mpos)
        const ChVector& GetPos()
        void SetPos_dt (const ChVector &mposdt)
        const ChVector& GetPos_dt()
        void SetPos_dtdt (const ChVector &mposdtdt)
        const ChVector& GetPos_dtdt()

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
        void SetForce(ChVector mf)
        ChVector& GetForce() 

    cdef cppclass ChNodeFEAxyzD(ChNodeFEAxyz):
        const ChVector& GetD()

    cdef cppclass ChNodeFEAxyzDD(ChNodeFEAxyzD):
        const ChVector& GetDD()

    cdef cppclass ChNodeFEAxyzrot(ChNodeFEAbase, ChBodyFrame):
        void SetForce(ChVector mf)
        ChVector& GetForce() 
        void SetTorque(ChVector mf)
        ChVector& GetTorque() 

    cdef cppclass ChNodeFEAxyzP(ChNodeFEAbase)

    cdef cppclass ChNodeFEAcurv(ChNodeFEAbase)


    # ------- ELEMENTS ------- #
    cdef cppclass ChElementBase:
        shared_ptr[ChNodeFEAbase] GetNodeN(int n)

    cdef cppclass ChElementGeneric(ChElementBase):
        ChElementGeneric()

    cdef cppclass ChElementBeam(ChElementGeneric):
        void EvaluateSectionDisplacement(const double eta, const ChMatrix &displ, ChVector &u_displ, ChVector &u_rotaz)
        void EvaluateSectionFrame(const double eta, const ChMatrix &displ, ChVector &u_displ, ChQuaternion &rot)
        void EvaluateSectionForceTorque(const double eta, const ChMatrix &displ, ChVector &Fforce, ChVector &Mtorque)
        void EvaluateSectionStrain(const double eta, const ChMatrix &displ, ChVector &StrainV)
        double GetMass()
        double GetRestLength()
        void SetRestLength(double ml)

    cdef cppclass ChElementCableANCF(ChElementBeam):
        void SetNodes(shared_ptr[ChNodeFEAxyzD] nodeA, shared_ptr[ChNodeFEAxyzD] nodeB)
        # void SetSection()
        void SetNodes(shared_ptr[ChNodeFEAxyzD] nodeA, shared_ptr[ChNodeFEAxyzD] nodeB)

    cdef cppclass ChElementBeamANCF:
        void SetNodes(shared_ptr[ChNodeFEAxyzDD] nodeA, shared_ptr[ChNodeFEAxyzDD] nodeB, shared_ptr[ChNodeFEAxyzDD] nodeC)
        void SetDimensions(double lenX, double beam_h, double beam_w)
        shared_ptr[ChNodeFEAxyzDD] GetNodeA()
        shared_ptr[ChNodeFEAxyzDD] GetNodeB()
        shared_ptr[ChNodeFEAxyzDD] GetNodeC()
        void setAlphaDamp(double a)
        double GetLengthX() const
        double GetMass()

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
        void SetAutomaticGravity(bool mg, int num_points=1)

    cdef cppclass ChMaterialSurface:
        ChMaterialSurface() except +

    cdef cppclass ChMaterialSurfaceSMC(ChMaterialSurface):
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

    cdef cppclass ChContactSurface:
        ChMesh* GetMesh()

    cdef cppclass ChContactSurfaceNodeCloud(ChContactSurface):
        ChContactSurfaceNodeCloud()
        void AddNode(shared_ptr[ChNodeFEAxyz] mnode, const double point_radius=0.001)
        void AddNode(shared_ptr[ChNodeFEAxyzrot] mnode, const double point_radius=0.001)
        void AddAllNodes(const double point_radius)

    cdef cppclass ChLinkPointFrame:
        int Initialize(shared_ptr[ChNodeFEAxyz] node, shared_ptr[ChBodyFrame] body, ChVector* pos)
        ChVector GetReactionOnNode()
        ChVector GetReactionOnBody()
        #virtual 
        #int Initialize(shared_ptr[ChNodeFEAxyz] node, shared_ptr[ChBodyFrame] body, ChVector[double] *pos=0)

    cdef cppclass ChLinkPointPoint:
        #virtual 
        ChLinkPointPoint()
        int Initialize(shared_ptr[ChNodeFEAxyz] anodeA, shared_ptr[ChNodeFEAxyz] anodeB)

