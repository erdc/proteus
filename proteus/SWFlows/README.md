# Proteus -- SWFlow app

The SWFlow app contains a number of tests for the the 2D nonlinear Shallow Water Equations, and the 2D modified Green-Naghdi equations. Both models are
numerically solved using continuous, linear finite elements as seen in
[Guermond, *et al* 2018](https://doi.org/10.1137/17M1156162) and [Guermond *et al* 2019].

The Shallow Water equations are a set of partial differential equations that form a hyperbolic system. They can be used to describe a body of water evolving under the action of gravity under the assumption that the deformations of the free surface are small compared to the water height.

The modified Green-Naghdi equations are a set of partial differential equations that form a
hyperbolic system and are an O(h) approximation to the traditional Green-Naghdi equations, where h is the mesh size. The Green-Naghdi equations are used to describe dispersive water waves.

The different tests demonstrate the respective models abilities to handle wet/dry states, reciprocate "real world" problems (see malpasset, colorado_river),
propagate solitary waves, etc.

## Getting Started

These instructions will help you run these tests on your local machine and help
you set up new tests for your specific needs.

### Prerequisites

The only prerequisite that you need is a working copy of Proteus. If you have yet to build Proteus, go to your Proteus directory and run

```
make develop
```

## Running the tests

Each test directory in `proteus/SWFlow/tests/` contains a python file named
`SWFlow.py`. This is the file that sets up the computational domain, initial conditions,
boundary conditions, and other problem specific options. Note that using the filename `SWFlow.py` is not necessary, one can use a different filename if wanted. Thus, in all following instructions, one can replace `SWFlow.py` with `some_filename.py`.

Each `SWFlow.py` file has its own set of default context options. These include choosing which model to use, final run time of the simulation, refinement level of mesh, etc. For example, in `SWFlow/tests/dam_over_bumps`, we see the following set of default options:

```
opts= Context.Options([
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time",30.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("refinement",4,"Level of refinement"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("reflecting_BCs",True,"Use reflecting BCs")
    ])

```
Here `sw_model` refers to which model is chosen: 0 for the Shallow Water Equations and 1 for the modified Green-Naghdi equations (or "Dispersive SWEs"). The rest of the options are straight forward.

To run a single test, go to a test directory of your choice. For this description, we will be using the directory `SWFlow/tests/dam_over_bumps` as an example. Without any modifications to the default context options, one can run the following command:

```
parun --SWEs -l1 -v SWFlow.py
```

which will execute the solver with the above default context options. The `-l1` flag controls the amount of output during the execution; you can choose any number between 1 and 10. The `-v` flag stands for verbose. (See the `proteus/parun.py` file for more details on the different flag options for `parun`.)

To use options other than the default options, one can run something like:

```
parun --SWEs -l1 -v SWFlow.py -C "sw_model=1 final_time=100"
```

which switches the model to the dispersive model and increases the final time to 100 seconds.

## Creating a test

Creating a new test is as simple as copying one of the existing directories,
and modifying appropriately to fit your needs.

Currently, most of the tests are set up to use a rectangular domain with triangle elements created by an in house function in Proteus. The SWFlow Apps can be used with meshes created by Triangle and Gmsh. The "real world" tests such as the Malpesset problem and Colorado River problem are using meshes created by the software SMS (https://www.aquaveo.com/software/sms-surface-water-modeling-system-introduction) in the ADH format. If you would like to use your own mesh, but aren't sure how to use it, please contact us and we can gladly help.

## Tests that aren't ready
elliptic_shoal
semicircular_shoal
colorado_river
south_padre_island
dunlap_dam
<!-- ## Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
``` -->

<!-- ## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds -->

## Authors

* **Dr. Manuel Quezada De Luna**
* **Dr. Chris Kees**
* **Eric Tovar(*)**

## Contributing

Please read [CONTRIBUTING.md](https://github.com/erdc/proteus/blob/master/CONTRIBUTING.md) for details on the Proteus code of conduct, and the process for submitting pull requests to us.

See also the list of Proteus [contributors](https://github.com/erdc/proteus/blob/master/CONTRIBUTORS.md) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

* J.-L. Guermond, M. Quezada de Luna, B. Popov, C. Kees, and M. Farthing. Well-balanced second-order finite element approximation of the shallow water equations with friction. SIAM Journal on Scientific Computing, 40(6):A3873– A3901, 2018.
* J.-L. Guermond, B. Popov, E. Tovar, C. Kees. Robust explicit relaxtion technique for solving the Green-Naghdi equations.
manuscript, 2019.

<!-- ## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc -->
