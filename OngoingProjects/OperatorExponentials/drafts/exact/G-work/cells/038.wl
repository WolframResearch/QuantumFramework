{a, ad} = FieldVariables[{1}];
g = BosonicMatrixElement[{m, m + 1}, a];
u = MatrixExp[-I lt {{0, g}, {g, 0}}];
claimed = {{Cos[lt Sqrt[m + 1]], -I BosonicMatrixElement[{m, m + 1}, a] Sin[lt Sqrt[m + 1]]/Sqrt[m + 1]},
   {-I BosonicMatrixElement[{m + 1, m}, ad] Sin[lt Sqrt[m + 1]]/Sqrt[m + 1], Cos[lt Sqrt[m + 1]]}};
{g, Simplify[ExpToTrig[u] - claimed, m >= 0]}
