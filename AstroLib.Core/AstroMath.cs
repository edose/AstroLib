using System.Diagnostics.CodeAnalysis;
using System;

namespace AstroLib.Core;

/// <summary>General math functions. Named maths plural not from any love of the British term,
/// but to distinguish this class from System.Math.</summary>
public static class AstroMath {

    /// <summary> True Euclidian modulo, always in range [ 0, |y| ), or null for y==0.</summary>
    /// <param name="x">numerator of modulo fraction.</param>
    /// <param name="y">denominator of modulo fraction. (Its sign is irrelevant for true Euclidian.)</param>
    /// <returns>True Euclidian modulo. Null if y==0. E.g. (7,3)=>1, (7,-3)=>1, (-7,3)=>2.</returns>
    public static double? EuclidianModulo(double x, double y) {
        if (y == 0) { return null; }
        return x - Math.Abs(y) * Math.Floor(x / Math.Abs(y));
    }
    
    /// <summary> Returns input value wrapped to within range [InputMin, InputMax).</summary>
    /// <param name="input">The value to be wrapped.</param>
    /// <param name="inputMin">Minimum of the range within which input shall be wrapped.</param>
    /// <param name="inputMax">Maximum of the range within which input shall be wrapped.
    /// Must be greater than inputMin, else exception thrown. </param>
    /// <returns>Input value as wrapped. </returns>
    public static double Wrap(double input, double inputMin, double inputMax) {
        if (inputMax - inputMin < 0) {
            throw new ArgumentException("Wrap() requires that inputMax > inputMin.");
        }
        if (input >= inputMin && input <= inputMax) 
            return input; // preserve exactly if within range.
        return EuclidianModulo(input - inputMin, inputMax - inputMin) + inputMin ?? 0.0;
    }

    /// <summary> Returns input value clamped to within range [InputMin, InputMax].
    /// Like System.Math.Clamp(), except allows null for each limit; each is ignored if null.
    /// Thus AstroMath.Clamp(x, null, null) always returns x.</summary>
    /// <param name="input">The value to be clamped.</param>
    /// <param name="inputMin">Minimum of the range within which input shall be clamped.</param>
    /// <param name="inputMax">Maximum of the range within which input shall be clamped.
    /// Must be greater than or equal to inputMin, else exception is thrown. </param>
    /// <returns>Input value as clamped </returns>
    public static double Clamp(double input, double? inputMin, double? inputMax) {
        var lowerLimit = inputMin ?? input;
        var upperLimit = inputMax ?? Math.Max(input, lowerLimit);
        return Math.Min(Math.Max(input, lowerLimit), upperLimit);
    }

    /// <summary> Returns true iff trial value is within valid range.
    /// If range's min or max is given as null, that limit does not apply.</summary>
    /// <param name="trialInput">The value to be tested.</param>
    /// <param name="inputMin">Minimum of the range against which input shall be tested.</param>
    /// <param name="inputMax">Maximum of the range against which input shall be tested.</param>
    /// <returns>True iff trial value is within valid range, for range limits not null;
    /// else false.</returns>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public static bool IsInRange(double trialInput, double? inputMin, double? inputMax) {
        return (AstroMath.Clamp(trialInput, inputMin, inputMax) == trialInput);
    }    
}