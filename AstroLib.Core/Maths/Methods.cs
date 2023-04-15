using static System.Math;

namespace AstroLib.Core.Maths;

/// <summary>This class simply holds static methods etc, thought of at random and
/// held here so as to not lose them.
/// Eventually all to be ported to the correct (future) class, somewhere in AstroLib.</summary>
public static class Methods {

    /// <summary> True Euclidian modulo, always in range [ 0, |y| ), or null for y==0.</summary>
    /// <param name="x">numerator of modulo fraction.</param>
    /// <param name="y">denominator of modulo fraction. (Its sign is irrelevant for true Euclidian.)</param>
    /// <returns>True Euclidian modulo. Null if y==0. E.g. (7,3)=>1, (7,-3)=>1, (-7,3)=>2.</returns>
    public static double? EuclidianModulo(double x, double y) {
        if (y == 0) { return null; }
        return x - Abs(y) * Floor(x / Abs(y));
    }
    
}