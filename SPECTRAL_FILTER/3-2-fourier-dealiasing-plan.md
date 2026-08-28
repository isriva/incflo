# 3/2 Fourier Dealiasing Plan

## Summary

Update `SpectralVelProductDecomp` to:

1. Forward-transform the original velocity with AMReX `FFT::R2C`.
2. Zero-pad/remap coefficients to a 3/2-sized spectral layout.
3. Inverse-transform once to a padded real-space velocity field.
4. Compute the outer product on the padded grid.
5. Forward-transform each product component on the padded grid.
6. Apply the existing `kmin`/`kmax` filter.
7. Truncate/remap coefficients to the original spectral layout.
8. Inverse-transform into the existing `vv_filter`.

## Implementation

- Reuse AMReX `FFT::R2C`, `forward`, `backward`, `getSpectralDataLayout`, and complex `FabArray` types.
- AMReX has no existing 3/2 dealiasing or spectral zero-padding/remapping utility. Its padding support is for OpenBC internals and is not applicable here.
- Add one small local helper for signed Fourier-mode remapping:
  - Initialize the destination spectrum to zero.
  - Copy modes using signed indices so negative frequencies move from the original upper range to the padded upper range.
  - Support both original-to-padded and padded-to-original transfers.
  - Preserve AMReX’s distributed spectral layouts.
- Add one hard-coded padding factor, defaulting to `1.5`.
- Use exact dimensions `3*N/2` in every active dimension and assert that each original dimension is even.
- Apply the padded FFT normalization when inverse-transforming padded velocity and when sending padded product coefficients to the original-grid inverse transform.
- Keep the function signature, output dimensions, component ordering, and current integer-index cutoff semantics unchanged.
- Support both 2D and 3D builds.
- Do not use `IncfloStructFact` for this path.

## Validation Constraint

Do not run tests, builds, formatters, or the application as part of this change.
