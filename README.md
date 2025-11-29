AlphaMatrixStuff: Converts fault current and voltages between phase domain and sequence domain. For vectorial addition, use Sum(Vtn) or Sum(Vabc)... for the majority of problems in EE45500, the vector sums should be zero as a sanity check.

ConvertingImpedance: Simple script to change base in per unit systems if a component does not match with the specified Vbase and/or Sbase.

NRaph: Load Flow project using the Newton Raphson method for a more efficient iterative process than Gauss Seidel method. Use Power Injection mismatch and Jacobian matrix to create updated guesses. This usually produces convergence in very low iteration counts. I plan on updating this script with MATPower to improve initializations and allow for handling of higher bus count system. 
