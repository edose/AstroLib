// ReSharper disable InconsistentNaming
namespace AstroLib.Observer; 

// TODO: Move this to new place?: AstroLib.Core, "Things" folder?

/// <summary> Represents a non-physical, abstract (numerically-defined) passband, e.g., SR for Sloan R.
///
/// To be precise, a passband has nothing to do with a filter, even if they share exactly the same name.
/// Often, in practice, a filter and passband are closely associated, e.g., flux data taken in
/// Sloan R filter is almost always transformed to magnitudes in standard Sloan R passband.
/// But this is not always true--for example, there is no consensus about which passband should be used
/// to express photometric data taken in a Clear filter--there is no "Clear" passband.
///
/// Transforms are used to convert flux data taken with a physical filter (piece of glass) to magnitudes
/// expressed in a standard (abstract) passband.</summary>
public class Passband {
    /// <summary> Identifier for this passband, e.g., "SR" or "V".</summary>
    public string ID { get; init; }
    
    /// <summary> CONSTRUCTOR. </summary>
    /// <param name="id">Identifier for this passband, e.g., "SR" or "V".</param>
    public Passband(string id) { ID = id; }
}
