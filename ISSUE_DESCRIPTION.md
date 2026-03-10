# 📝 Add Comments & Documentation to Key Algorithms

## Issue
Several critical algorithms in the mapping module lack sufficient inline documentation, hindering developer onboarding and code maintenance.

## Problem
- **PartitionOfUnityMapping**: Complex clustering and weight normalization not explained
- **RadialBasisFctSolver**: Matrix assembly and decomposition strategies underdocumented  
- **Dead Axis Handling**: Numerical conditioning for 2D/3D data lacks clarity

## Solution
Add comprehensive inline documentation explaining:

### `PartitionOfUnityMapping.hpp`
- PUM algorithm overview and workflow
- Shepard's weight normalization method
- Conservative vs. consistent mapping strategies

### `RadialBasisFctSolver.hpp`
- RBF matrix assembly (C and A matrices)
- Decomposition methods (Cholesky vs. QR)
- Polynomial basis function handling
- Dead axis detection and conditioning

## Implementation Details
- Documentation added as code comments (no external files)
- Includes mathematical formulas and algorithm steps
- Explains numerical stability considerations

## Files Changed
- `src/mapping/PartitionOfUnityMapping.hpp` (+65 lines)
- `src/mapping/RadialBasisFctSolver.hpp` (+132 lines)

## Benefits
✅ Improves code readability  
✅ Simplifies debugging  
✅ Supports code reviews  
✅ Reduces onboarding time  

## Type
- Documentation/Code Quality
- Difficulty: Beginner-friendly
- Effort: ~3 hours
