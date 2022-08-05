# Don't pollute the global namespace
def _do_speedup():
    def _get_version_number():
        """
        version() returns a string (e.g.)
        'SageMath version 9.6, Release Date: 2022-05-15'

        The below parses the string and returns the version as a float
        float('9.6')
        """
        return float(version().split(' ')[2][:-1])

    # And because why not
    proof.all(False)

    # Lorenz Panny has fixed both of the below monkey patches with the tickets:
    # - https://trac.sagemath.org/ticket/34281 (Caching of the finite fields)
    # - https://trac.sagemath.org/ticket/34284 (Dimension of hyperelliptic curve)
    #
    # We should check the version of sage and if >= 9.7 skip the below patches
    if _get_version_number() < 9.7:
        # Since this type gets created before we could ever hope to monkey patch the 
        # `sage.categories.fields.Fields.ParentMethods`
        # method, we'll patch it on the relevant type instead.
        # We'll patch a few different types to make sure we get the relevant things (large and small prime, extension and no extension)
        p = 2^127 - 1 # Arbitrary large prime
        to_patch = [GF(3), GF(3^2), GF(p), GF(p^2)]
        for x in to_patch:
            type(x).vector_space = sage.misc.cachefunc.cached_method(type(x).vector_space)

        # An alternative would be to replace the bytecode in 
        # `sage.categories.fields.Fields.ParentMethods.vector_space`
        # as all types share the same method, by identity
        # Something to be explored later, perhaps :)
        
        # No use calculating the dimension of HyperElliptic every single time
        from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
        AlgebraicScheme_subscheme_projective.dimension = lambda self: 1
    

_do_speedup()
