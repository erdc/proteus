# How to run it?

branch: yy/use_particle_velocity_in_c
commit: 0f539f22495f45033930e5127f1cf47a2f1994ec



# Right now rho=1.5 begins to work
# Instability when two particles get close to each other
This is because `dT_Chrono=ct.dt_fixed/1000.0`.
Use `dT_Chrono=ct.dt_fixed/100.0` to get good result.
# Turn off the collison in Chrono
auto msphereBody = std::make_shared<ChBodyEasySphere>(ball_radius[i],        // radius size
                                                                  m_particles_density, // density
                                                                  false,                // collide enable?
                                                                  false);              // visualization?
# Use repulsive force model in Glowinski's paper
# VBDF is totally wrong.
#######################################################################################################################################
05/31/2018
#1. unstructured mesh is not good
#2. turn on supg is not good
#3. there is no change to use particle_penalty_constant=1e3 instead of 1e6
#4. use `PSTAB=0.0` has a little improvement
#5. Repulsive force model (19-20) in Glowinski et. al. paper is similar to 
chrono's result with default setting.
#6. No big improvement using `particle_epsFact=2.0` instead of `1.0`
#7. Using penalty method to get pressure in .h gives blowup of the pressure.
#7. Use entropy viscosity, ` self.ARTIFICIAL_VISCOSITY = 2`, it is not good; sudden change;
 
#######################################################################################################################################
# It is because of `self.model.q[('velocityStar', 0)]` the particle changes the motion suddently.
#######################################################################################################################################
05/30/2018
# 1. The behavior of `p_sharp` with `v_star=vn` is not good. So use `p_sharp=pn` as follows in `pressure_p.py`
```
    def postStep(self, t, firstStep=False):
        """
        Give the TC object an opportunity to modify itself before the time step.
        """
 
        self.model.q_p_sharp[:] = self.model.q[('u', 0)]
        self.model.ebqe_p_sharp[:] = self.model.ebqe[('u', 0)]
        self.model.q_grad_p_sharp[:] = self.model.q[('grad(u)', 0)]
        self.model.ebqe_grad_p_sharp[:] = self.model.ebqe[('grad(u)', 0)]
         
        copyInstructions = {}
        return copyInstructions
```

# 2. Turn on rotation term in pressure
There is no difference in this case.

#3. Is this an eError: `c[('f', 0)][..., i] = np.min(rho * nu) * velocity[..., i]` in `Pres.h`
(There is no `rho` here )?
No. `mu=rho*nu`

#4. The update of velocity at quadrature point should be after the computation of pressure
since the rotation form should use `\tilde{u}`

