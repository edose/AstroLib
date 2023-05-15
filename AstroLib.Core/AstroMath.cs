using System.Diagnostics.CodeAnalysis;

namespace AstroLib.Core;

/// <summary>Basic math functions frequently useful for astronomy.</summary>
public static class AstroMath {
    /// <summary> True Euclidian modulo, always non-negative and always in range [ 0, |y| )
    /// (Note: excludes the exact top of range in favor of the floor value). Value is null for y==0.
    /// The best modulo for wrapping values (e.g., longitudes) than the floor or other modulo
    /// functions, as it is unambiguous and well-behaved near the range min and max values,
    /// and whenever y is negative, and as it can never yield negative values.</summary>
    /// <param name="x">numerator of modulo fraction.</param>
    /// <param name="y">denominator of modulo fraction. (Its sign is irrelevant for true Euclidian.)</param>
    /// <returns>True Euclidian modulo. Value is null if y==0.\n
    /// &emsp;`AstroMath.EuclidianModulo(7, 3)` => 1\n
    /// &emsp;`AstroMath.EuclidianModulo(7, -3)` => 1\n
    /// &emsp;`AstroMath.EuclidianModulo(-7, 3)` => 2\n
    /// &emsp;`AstroMath.EuclidianModulo(x, 0)` => null, for any x.\n</returns>
    public static double? EuclidianModulo(double x, double y) {
        if (y == 0) { return null; }
        return x - Math.Abs(y) * Math.Floor(x / Math.Abs(y));
    }
    
    /// <summary> Returns input value as wrapped to within range [ InputMin, InputMax ), *excluding*
    /// the range maximum. Especially useful for managing longitudes.</summary>
    /// <param name="input">The value to be wrapped.</param>
    /// <param name="rangeMin">Minimum of the range within which input shall be wrapped.</param>
    /// <param name="rangeMax">Maximum of the range within which input shall be wrapped.
    /// Must be greater than rangeMin, else ArgumentException thrown. </param>
    /// <returns>Input value as wrapped.\n
    /// &emsp;`AstroMath.Wrap( -7, 0, 360)` => 353\n
    /// &emsp;`AstroMath.Wrap(  0, 0, 360)` =>   0\n
    /// &emsp;`AstroMath.Wrap( 77, 0, 360)` =>  77\n
    /// &emsp;`AstroMath.Wrap(360, 0, 360)` =>  77\n
    /// &emsp;`AstroMath.Wrap(367, 0, 360)` =>   7\n
    /// &emsp;`AstroMath.Wrap(-455, -180, 180)` =>  -95\n
    /// &emsp;`AstroMath.Wrap(-180, -180, 180)` => -180\n
    /// &emsp;`AstroMath.Wrap(   7, -180, 180)` =>    7\n
    /// &emsp;`AstroMath.Wrap( 180, -180, 180)` => -180 (note: max is excluded)\n
    /// &emsp;`AstroMath.Wrap( 455, -180, 180)` =>   95\n
    /// &emsp;`AstroMath.Wrap( 1.333, -180, 180)` =>   1.333\n
    ///Note: as in the last example above, if input value is already within range,
    /// it is returned directly so that floating-point equality tests remain valid.</returns>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public static double Wrap(double input, double rangeMin, double rangeMax) {
        if (rangeMax - rangeMin < 0) {
            throw new ArgumentException("Wrap() requires that inputMax > inputMin.");
        }
        if (input >= rangeMin && input < rangeMax) 
            return input; // preserve exactly if within range, inputMax excluded.
        if (input == rangeMax)
            return rangeMin;  // range excludes exact inputMax.
        return EuclidianModulo(input - rangeMin, rangeMax - rangeMin) + rangeMin ?? 0.0;
    }

    /// <summary> Returns input value clamped to within range [ InputMin, InputMax ],
    /// inclusive of both range limits. Useful for handling latitudes, declinations, and polar angles.
    /// Like native System.Math.Clamp(), except AstroMath.Clamp allows null for each limit;
    /// each is ignored if null. AstroMath.Clamp(x, null, null) always returns x.</summary>
    /// <param name="input">The value to be clamped.</param>
    /// <param name="rangeMin">Minimum of the range within which input shall be clamped.</param>
    /// <param name="rangeMax">Maximum of the range within which input shall be clamped.
    /// Must be greater than or equal to rangeMin, else exception is thrown. </param>
    /// <returns>Input value clamped to given range.\n
    /// &emsp;`AstroMath.Clamp(-100, -23, 31)` => -23\n
    /// &emsp;`AstroMath.Clamp(-23, -23, 31)` => -23\n
    /// &emsp;`AstroMath.Clamp(1, -23, 31)` => 1\n
    /// &emsp;`AstroMath.Clamp(23, -23, 31)` => 31\n
    /// &emsp;`AstroMath.Clamp(100, -23, 31)` => 31\n
    /// &emsp;`AstroMath.Clamp(x, null, null)` => x, for any x.\n
    /// &emsp;`AstroMath.Clamp(1.333, -23, 31)` => 1.333\n
    /// Note: as in the last example above, if input value is already within range,
    /// it is returned directly so that floating-point equality tests remain valid.</returns>
    public static double Clamp(double input, double? rangeMin, double? rangeMax) {
        var lowerLimit = rangeMin ?? input;
        var upperLimit = rangeMax ?? Math.Max(input, lowerLimit);
        return Math.Min(Math.Max(input, lowerLimit), upperLimit);
    }

    /// <summary> Returns true iff trial value is within range [ rangeMin, rangeMax], inclusive of limits.
    /// If range's min or max is given as null, that limit does not apply.</summary>
    /// <param name="trialInput">The value to be tested.</param>
    /// <param name="rangeMin">Minimum of the range against which input shall be tested.</param>
    /// <param name="rangeMax">Maximum of the range against which input shall be tested.</param>
    /// <returns>True iff trial value is within given range, for range limits not null; else false.\n
    /// &emsp;`AstroMath.IsInRange(-2, 5, 15)` => false\n
    /// &emsp;`AstroMath.IsInRange(5, 5, 15)` => true\n
    /// &emsp;`AstroMath.IsInRange(6, 5, 15)` => true\n
    /// &emsp;`AstroMath.IsInRange(15, 5, 15)` => true\n
    /// &emsp;`AstroMath.IsInRange(22, 5, 15)` => false\n
    /// &emsp;`AstroMath.IsInRange(x, null, null)` => true, for any x.</returns>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public static bool IsInRange(double trialInput, double? rangeMin, double? rangeMax) {
        return (AstroMath.Clamp(trialInput, rangeMin, rangeMax) == trialInput);
    }

    /// <summary> Return new 2-dimensional array of doubles from 1-dimensional array
    /// (without changing memory layout or total number of elements).</summary>
    /// <param name="array1d">One-dimensional array supplying values</param>
    /// <param name="size1">length of 2-D array in first index</param>
    /// <param name="size2">length of 2-D array in second index</param>
    /// <returns>2-Dimensional array of doubles</returns>
    /// <exception cref="ArgumentException">Thrown if size1*size2 != array1d.Length</exception>
    public static double[,] Reshape1dTo2dArray(double[] array1d, int size1, int size2) {
        if (array1d.Length != size1 * size2) {
            throw new ArgumentException("array1d.Length must equal size1 * size2 but does not.");
        }
        var array2d = new double[size1, size2];
        for (int i = 0; i < size1 * size2; ++i) {
            array2d[i / size2, i % size2] = array1d[i];
        }
        return array2d;
    }
    
    
}