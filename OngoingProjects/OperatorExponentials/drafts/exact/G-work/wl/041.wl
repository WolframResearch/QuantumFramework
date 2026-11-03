sx = {{0, 1}, {1, 0}}; u = MatrixExp[-I lam t m sx];
{Simplify[ExpToTrig[u - (Cos[lam t m] IdentityMatrix[2] - I Sin[lam t m] sx)]],
 Simplify[ExpToTrig[u - (Outer[Times, {1, 1}, {1, 1}]/2 Exp[-I lam t m] + Outer[Times, {1, -1}, {1, -1}]/2 Exp[I lam t m])]]}
