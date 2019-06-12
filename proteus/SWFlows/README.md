# Proteus -- SWFlows app

The SWFlows app contains a number of tests for the the 2D nonlinear Shallow Water Equations, and the 2D modified Green-Naghdi equations. Both models are
numerically solved using continuous, linear finite elements as seen in
[cite friction paper] and [cite mGN paper].

Some introduction stuff on the SW equations.

The modified Green-Naghdi equations are a set of a equations that form a
hyperbolic system and are an O(h) approximation to the traditional Green-Naghdi
equations, where h is the mesh size.  

The different tests demonstrate the respective models abilities to handle wet/dry states, reciprocate "real" world problems (see malpasset, colorado_river),
propagate solitary waves, etc etc etc.

## Getting Started

These instructions will help you run these tests on your local machine and help
you set up new tests for your specific needs.

### Prerequisites

Each test directory in `/SWFlows/tests/` contains a python file named
SWFlow.py. This is the file that sets up a domain, initial conditions and
boundary conditions for a specific problem. This is the only file that you
need to modify to run these tests.

Each SWFlow.py file has its own set of options to choose which model to run,
final run time of the simulation, refinement level of mesh, etc. For example,
in `SWFlows/tests/dam_over_bumps`, we see the following set of default options:

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


## Running the tests

To

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Manuel Quezada De Luna**
* **Chris Kees**
* **Eric Tovar**

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

* J.-L. Guermond, M. Quezada de Luna, B. Popov, C. Kees, and M. Farthing. Well-balanced second-order finite element approximation of the shallow water equations with friction. SIAM Journal on Scientific Computing, 40(6):A3873– A3901, 2018.
* J.-L. Guermond,, B. Popov, E. Tovar, C. Kees. Robust explicit relaxtion technique for solving the Green-Naghdi equations.
manuscript, 2019

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
