# Charged-current DIS and SIDIS

## Pitfalls

### Mysterious template substitution failure errors

Check that you are calling the correct function. Many internal functions take a `auto &` parameter and assume that you have passed in a correct type. For example, `void SIDIS::lepton_pair_xy` expects an argument of type `LHAInterface` for the parameter `pdf` while `void SIDIS::lepton_pair_xy_errors` expects an `LHASetInterface`. Mixing these will result in a template substitution error, because the function is trying to access some property of the `pdf` object that doesn't exist.

Remember to always check from the error message which line ultimately triggers the substitution failure.

### No member named '...' in some internal type

See Mysterious template substitution failure errors.
