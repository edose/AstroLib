using System.Diagnostics.CodeAnalysis;

namespace AstroLib.Core.Geometry; 

/// <summary>Represents a plane angle.</summary>
public class Angle : IEquatable<Angle> {
    public double Radians { get; init; }
    public double Degrees    => Radians * Constants.DegreesPerRadian;
    public double Arcminutes => Radians * Constants.ArcminutesPerRadian;
    public double Arcseconds => Radians * Constants.ArcsecondsPerRadian;

    /// <summary>Private default constructor, to prevent user construction without input data.</summary>
    private protected Angle() {}

    /// <summary> Constructor from a double representing radians. </summary>
    /// <param name="radians">Magnitude of angle, in radians.</param>
    public Angle(double radians) { Radians = radians; }

    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <returns>Angle instance.</returns>
    public static Angle FromDegrees(double degrees) {
        return new Angle(degrees / Constants.DegreesPerRadian);
    }

    /// <summary> Returns true iff other Angle exactly equals this one in magnitude.</summary> 
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public bool Equals(Angle? other) {
        if (other == null) { return false; }
        return (other.Radians == Radians);
    }
    
    public override string ToString() {
        return $"Angle: {Radians} radians = {Degrees} degrees.";
    }
}

/// <summary>Represents a Latitude.</summary>
public sealed class Latitude : Angle {
    public readonly double ClampMin = -Math.PI / 2.0;
    public readonly double ClampMax = +Math.PI / 2.0;

    /// <summary>Private default constructor, to prevent construction without input data.</summary>
    private Latitude() {}

    /// <summary>Constructor from a double representing radians.</summary>
    /// <param name="radians">The magnitude of this Latitude, limited to [-pi/2, pi/2].</param>
    /// <param name="clampToRange">If true, clamp values outside the validity range [-pi/2, pi/2]
    /// to that range. If false, values outside the validity range cause an exception.</param>
    /// <exception cref="ArgumentException">Thrown if input value is outside the validity range and
    /// user has disallowed clamping (via clampToRange=false).</exception>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public Latitude(double radians, bool clampToRange = false) : base(radians) {
        var clampedRadians = AstroMath.Clamp(radians, ClampMin, ClampMax);
        if (clampedRadians != radians && !clampToRange) {
            throw new ArgumentException("Latitude value (radians) is out of range " +
                                        "and clamping to valid range is disabled.");
        }
        Radians = clampedRadians;
    }

    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <param name="clampToRange"></param>
    /// <returns>Latitude instance.</returns>
    public new static Latitude FromDegrees(double degrees, bool clampToRange = false) {
        return new Latitude(degrees / Constants.DegreesPerRadian, clampToRange);
    }

    /// <summary> Returns new PolarAngle object that is equivalent to this Latitude object.</summary>
    public PolarAngle ToPolarAngle() {
        return new PolarAngle(Math.PI / 2.0 - Radians, true);
    }

    public override string ToString() {
        return $"Latitude: {Radians} radians = {Degrees} degrees.";
    }
}

/// <summary>Represents a Polar Angle, a spherical coordinate's angle with the north pole.
/// Preferred by some mathematicians over Latitude (which is preferred in astronomy).</summary>
public sealed class PolarAngle : Angle {
    public readonly double ClampMin = 0;
    public readonly double ClampMax = Math.PI;
    
    /// <summary>Private default constructor, to prevent construction without input data.</summary>
    private PolarAngle() {}
    
    /// <summary>Constructor from a double representing radians.</summary>
    /// <param name="radians">The magnitude of this Polar Angle, limited to [0, pi].</param>
    /// <param name="clampToRange">If true, clamp values outside the validity range [0, pi]
    /// to that range. If false, values outside the validity range cause an exception.</param>
    /// <exception cref="ArgumentException">Thrown if input value is outside the validity range and
    /// user has disallowed clamping (via clampToRange=false).</exception>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public PolarAngle(double radians, bool clampToRange=false): base(radians){
        var clampedRadians = AstroMath.Clamp(radians, ClampMin, ClampMax);
        if (clampedRadians != radians && !clampToRange) {
            throw new ArgumentException("PolarAngle value (radians) is out of range " +
                                        "and clamping to valid range is disabled.");
        }
        Radians = clampedRadians;
    }
    
    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <param name="clampToRange"></param>
    /// <returns>PolarAngle instance.</returns>
    public new static PolarAngle FromDegrees(double degrees, bool clampToRange = false) {
        return new PolarAngle(degrees / Constants.DegreesPerRadian, clampToRange);
    }
    
    /// <summary> Returns new Latitude object that is equivalent to this PolarAngle object.</summary>
    public Latitude ToLatitude() {
        return new Latitude(Math.PI / 2.0 - Radians, true);
    }
    
    public override string ToString() {
        return $"PolarAngle: {Radians} radians = {Degrees} degrees.";
    }
}
    
/// <summary>Represents a Longitude.</summary>
public sealed class Longitude : Angle {
    public readonly double WrapMin = -Math.PI;
    public readonly double WrapMax = +Math.PI;
    
    /// <summary>Private default constructor, to prevent construction without input data.</summary>
    private Longitude() {}
    
    /// <summary>Constructor from a double representing radians.</summary>
    /// <param name="radians">The magnitude of this Longitude, limited to [-pi, +pi).</param>
    /// <param name="wrapToRange">If true, wrap values outside the validity range [-pi, +pi) to inside
    /// that range. If false, values outside the validity range cause an ArgumentException.</param>
    /// <exception cref="ArgumentException">Thrown if input value is outside the validity range and
    /// user has disallowed clamping (via clampToRange=false).</exception>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public Longitude(double radians, bool wrapToRange=false) : base(radians) {
        var wrappedRadians = AstroMath.Wrap(radians, WrapMin, WrapMax);
        if (wrappedRadians != radians && !wrapToRange) {
            throw new ArgumentException("Longitude value (radians) is out of range " +
                                        "and wrapping into valid range is disabled.");
        }
        Radians = wrappedRadians;
    }
    
    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <param name="wrapToRange"></param>
    /// <returns>Longitude instance.</returns>
    public new static Longitude FromDegrees(double degrees, bool wrapToRange = false) {
        return new Longitude(degrees / Constants.DegreesPerRadian, wrapToRange);
    }
    
    public override string ToString() {
        return $"Longitude: {Radians} radians = {Degrees} degrees.";
    }
}

/// <summary>Represents a Declination, the sky analog to terrestrial Latitude, and generally
/// expressed in degrees within [-90, +90].</summary>
public sealed class Declination : Angle {
    public readonly double WrapMin = -Math.PI / 2.0;
    public readonly double WrapMax = +Math.PI / 2.0;

    /// <summary>Private default constructor, to prevent construction without input data.</summary>
    private Declination() {}

    /// <summary>Constructor from a double representing radians.</summary>
    /// <param name="radians">The magnitude of this Declination in radians,
    /// limited to [-pi/2, pi/2].</param>
    /// <param name="clampToRange">If true, clamp values outside the validity range [-pi/2, pi/2]
    /// to that range. If false, values outside the validity range cause an exception.</param>
    /// <exception cref="ArgumentException">Thrown if input value is outside the validity range and
    /// user has disallowed clamping (via clampToRange=false).</exception>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public Declination(double radians, bool clampToRange = false) : base(radians) {
        var clampedRadians = AstroMath.Clamp(radians, WrapMin, WrapMax);
        if (clampedRadians != radians && !clampToRange) {
            throw new ArgumentException("Declination value (radians) is out of range " +
                                        "and clamping to valid range is disabled.");
        }
        Radians = clampedRadians;
    }
    
    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <param name="clampToRange"></param>
    /// <returns>Declination instance.</returns>
    public new static Declination FromDegrees(double degrees, bool clampToRange = false) {
        return new Declination(degrees / Constants.DegreesPerRadian, clampToRange);
    }
    
    /// <summary> Returns new PolarAngle object that is equivalent to this Declination object.</summary>
    public PolarAngle ToPolarAngle() {
        return new PolarAngle(Math.PI / 2.0 - Radians, true);
    }
    
    public override string ToString() {
        return $"Declination: {Radians} radians = {Degrees} degrees.";
    }
}

// TODO: Longitude and RightAscension need to *wrap* input angles, not to clamp them.

/// <summary>Represents a RightAscension, the sky analog to terrestrial Longitude, and generally
/// expressed in hours within [0, 24).</summary>
public sealed class RightAscension : Angle {
    public readonly double WrapMin = 0;
    public readonly double WrapMax = 2 * Math.PI;
    public double Hours => Degrees / 15.0;

    /// <summary>Private default constructor, to prevent construction without input data.</summary>
    private RightAscension() {}

    /// <summary>Constructor from a double representing radians. By contrast to other </summary>
    /// <param name="radians">The magnitude of this Longitude in radians, limited to [0, 2*pi).</param>
    /// <param name="wrapToRange">If true, wrap values outside the validity range [0, 2*pi)
    /// to that range by adding or subtracting whatever required multiples of 2*pi.
    /// If false, values outside the validity range cause an exception.</param>
    /// <exception cref="ArgumentException">Thrown if input value is outside the validity range and
    /// user has disallowed clamping (via clampToRange=false).</exception>
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public RightAscension(double radians, bool wrapToRange=false) : base(radians) {
        var wrappedRadians = AstroMath.Wrap(radians, WrapMin, WrapMax);
        if (wrappedRadians != radians && !wrapToRange) {
            throw new ArgumentException("RightAscension value (radians) is out of range " +
                                        "and wrapping into valid range is disabled.");
        }
        Radians = wrappedRadians;
    }
    
    /// <summary> Factory method using degrees. </summary>
    /// <param name="degrees">Magnitude of angle, in degrees. </param>
    /// <param name="wrapToRange"></param>
    /// <returns>RightAscension instance.</returns>
    public new static RightAscension FromDegrees(double degrees, bool wrapToRange = false) {
        return new RightAscension(degrees / Constants.DegreesPerRadian, wrapToRange);
    }
    
    public override string ToString() {
        return $"RightAscension: {Radians} radians = {Degrees} degrees.";
    }
}
