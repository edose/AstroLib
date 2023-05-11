// ReSharper disable ClassNeverInstantiated.Global

namespace AstroLib.Observer;

// TODO: Move this to new place?: AstroLib.Core, "Things" folder?

/// <summary> Represents a transform, that is, an adjustment in magnitudes from an instrumental magnitude
/// taken through a specific filter to a magnitude in a specific passband. Requires a color index, which is
/// defined here as ColorPassband1 - ColorPassband2.</summary>
public class Transform {
    public Filter Filter;
    public Passband Passband;
    public Passband ColorPassband1;
    public Passband ColorPassband2;
    public double TransformValue; 
        
    /// <summary> CONSTRUCTOR for the normal case. </summary>
    /// <param name="filter"></param>
    /// <param name="passband"></param>
    /// <param name="colorPassband1"></param>
    /// <param name="colorPassband2"></param>
    /// <param name="transformValue"></param>
    public Transform(Filter filter, Passband passband, Passband colorPassband1, Passband colorPassband2,
        double transformValue) {
        Filter = filter;
        Passband = passband;
        ColorPassband1 = colorPassband1;
        ColorPassband2 = colorPassband2;
        TransformValue = transformValue;
    }
}