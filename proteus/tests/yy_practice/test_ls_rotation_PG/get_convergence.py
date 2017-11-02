import shelve
from contextlib import closing
import numpy
import sys


def get_convergence_rate(datafile):
    h = []
    error = {}
    level_convergence = {}
    convergence = {}

    nl = 0

    with closing(shelve.open(datafile)) as data:
        nl = len(data['mesh'])

        for l in range(nl):
            h += data['mesh'][l]['h']

            error_dict_l = data['errorData'][0][l]
            for error_type in error_dict_l:
                if error_type[:5] == 'error':
                    try:
                        error[error_type] += error_dict_l[error_type]
                    except KeyError:
                        error[error_type] = error_dict_l[error_type]
                        level_convergence[error_type] = []
                        convergence[error_type] = []
    #     print data

    # print h
    # print error

    for l in range(1, nl):
        for error_type in error:
            level_convergence[error_type].append(
                numpy.log(error[error_type][l] / error[error_type][l - 1]) / numpy.log(h[l] / h[l - 1]))

    # print level_convergence
    if nl > 1:
        for error_type in level_convergence:
            convergence[error_type] = sum(
                level_convergence[error_type]) / len(level_convergence[error_type])
    return error, convergence


if __name__ == '__main__':
    import rotation2D
    from proteus import Context, Profiling, Comm
    comm = Comm.init()
    Profiling.procID = comm.rank()

    if len(sys.argv) > 1:
        print get_convergence_rate(sys.argv[1])
    else:
        ctx = Context.get()  # To get filename
        print get_convergence_rate(ctx.datafile)
