
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)

cdef extern from "ChMoorings.h":
    cdef cppclass ChMesh:
        void SetAutomaticGravity(bool mg, int num_points=1)
    cdef cppclass ChMaterialSurfaceBase
    cdef cppclass ChMaterialSurfaceDEM:
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
    cdef cppclass ChContactSurfaceNodeCloud:
        ChContactSurfaceNodeCloud()
        void AddNode(shared_ptr[ChNodeFEAxyz] mnode, const double point_radius=0.001)
        void AddNode(shared_ptr[ChNodeFEAxyzrot] mnode, const double point_radius=0.001)
        void AddAllNodes(const double point_radius)
    cdef cppclass ChNodeFEAxyz
    cdef cppclass ChNodeFEAxyzrot
    cdef cppclass ChNodeFEAxyzDD:
        const ChVector& GetPos()
        const ChVector& GetD()
        const ChVector& GetDD()
        void SetPos (const ChVector &mpos)
        const ChVector& GetPos_dt()
        void SetPos_dt (const ChVector &mposdt)
        const ChVector& GetPos_dtdt()
        void SetPos_dtdt (const ChVector &mposdtdt)
        void SetFixed(bool mev)
    cdef cppclass ChElementBeamANCF:
        void SetNodes(shared_ptr[ChNodeFEAxyzDD] nodeA, shared_ptr[ChNodeFEAxyzDD] nodeB, shared_ptr[ChNodeFEAxyzDD] nodeC)
        void SetDimensions(double lenX, double beam_h, double beam_w)
        shared_ptr[ChNodeFEAxyzDD] GetNodeA()
        shared_ptr[ChNodeFEAxyzDD] GetNodeB()
        shared_ptr[ChNodeFEAxyzDD] GetNodeC()
        void setAlphaDamp(double a)
        double GetLengthX() const
    cdef cppclass ChBodyFrame
    cdef cppclass ChLinkPointFrame:
        ChVector GetReactionOnNode()
        #virtual 
        #int Initialize(shared_ptr[ChNodeFEAxyz] node, shared_ptr[ChBodyFrame] body, ChVector[double] *pos=0)
    cdef cppclass ChLinkPointPoint:
        #virtual 
        ChLinkPointPoint()
        int Initialize(shared_ptr[ChNodeFEAxyz] anodeA, shared_ptr[ChNodeFEAxyz] anodeB)
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
    cdef cppclass ChMatrix33[double]:
        ChVector Get_A_Xaxis()
        ChVector Get_A_Yaxis()
        ChVector Get_A_Zaxis()
    cdef cppclass ChSystem:
        void Add(shared_ptr[ChPhysicsItem] newitem)
    cdef cppclass ChSystemDEM:
        void Add(shared_ptr[ChPhysicsItem] newitem)
        double GetStep()
        double GetChTime()
        void SetupInitial()
    cdef cppclass ChPhysicsItem:
        ChPhysicsItem()
    cdef cppclass ChBody(ChPhysicsItem):
        ChBody() except +
        void SetPos(const ChVector[double]& mpos) except +
        void SetRot(ChQuaternion &rot) except +
        void SetBodyFixed(bool state) except +
        void SetMaterialSurface(const shared_ptr[ChMaterialSurfaceBase] &mnewsurf) except +
        void SetMass(double newmass)
    cdef cppclass ChBodyEasyBox(ChBody):
        ChBodyEasyBox(double Xsize, double Ysize, double Zsize, double mdensity, bool collide=False, bool visual_asset=True)
