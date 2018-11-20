Navier-Stokes
*************


Description
===========

There are currently 3 implementations of Navier-Stokes equations in proteus:

* Two-phase flow (air/water)
* Three-phase flow (air/water/sediment)

Two-Phase
=========
The two-phase implementation of Navier-Stokes.
Files corresponding to model:
`proteus.mprans.RANS2P.py`,

`proteus.mprans.RANS2P.h`,
`proteus.mprans.cRANS2P.pyx`

`proteus.mprans.RANS2P2D.h`,
`proteus.mprans.cRANS2P2D.pyx`


Three-Phase
===========
The three-phase implementation of Navier-Stokes.
Files corresponding to model:
`proteus.mprans.RANS3P.py`,
`proteus.mprans.RANS3PSed.py`,

`proteus.mprans.RANS3PF.h`,
`proteus.mprans.RANS3PF2D.h`,
`proteus.mprans.cRANS3PF.pyx`

`proteus.mprans.RANS3PSed.h`,
`proteus.mprans.RANS3PSed2D.h`,
`proteus.mprans.cRANS3PSed.pyx`


Dealing with Moving Domains
===========================

When dealing with moving domains, the option `movingDomain` must be set to `True`.


Moving (ALE) Mesh
-----------------

Immersed Boundaries
-------------------

`proteus.mprans.RANS2P.py`,

`proteus.mprans.RANS2P.h`,
`proteus.mprans.cRANS2P.pyx`
