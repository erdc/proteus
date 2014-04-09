def test_povich():
    """test_povich  tests the l_2 Norm of a vector

    Verifies that the proteus.LinearAlgebraTools.Vec constructor
    correctly creates one-dimensional arrays of the given length of
    type double precision and with entries set to zero for several
    trials.
    """
    from proteus.LinearAlgebraTools import l2Norm 
    for n in [1, 10, 100, 1000]:
        # Generate a vector of n zeros
        x = np.zeros(n, )
        # Vector of length n
        eq(x.size, n)
        # Norm of a Vector of Zeros is Zero
        eq(l2Norm(x)<= 1.0e-39, True)
