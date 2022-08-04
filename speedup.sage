# Don't pollute the global namespace
def _do_speedup():
    # And because why not
    proof.all(False)

    # Since this type gets created before we could ever hope to monkey patch the 
    # `sage.categories.fields.Fields.ParentMethods`
    # method, we'll patch it on the relevant type instead.
    # We'll patch a few different types to make sure we get the relevant things (large and small prime, extension and no extension)
    p = 2^127 - 1 # Arbitrary large prime
    to_patch = [GF(3), GF(3^2), GF(p), GF(p^2)]
    for x in to_patch:
        type(x).vector_space = sage.misc.cachefunc.cached_method(type(x).vector_space)



    sage.schemes.hyperelliptic_curves.hyperelliptic_generic.Hyperelliptic_generic = lambda self: 1

    # An alternative would be to replace the bytecode in 
    # `sage.categories.fields.Fields.ParentMethods.vector_space`
    # as all types share the same method, by identity
    # Something to be explored later, perhaps :)
    

_do_speedup()
