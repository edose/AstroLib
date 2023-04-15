namespace AstroLib.Core.Geometry; 

public static class Spherical {
    /// <summary>Returns list of spherical coordinates evenly spaced on a sphere.</summary>
    /// <param name="pointCount">Number of points wanted.</param>
    /// <returns>List of SphCoord objects, each with null radius (direction-only).
    ///     Note: polarAngle is not Latitude or Declination.</returns>
    public static List<SphCoord> MakeGoldenSpiral(int pointCount) {
        var points = new List<SphCoord>();
        for (var i = 0; i < pointCount; i++) {
            var index = i + 0.5;
            var azimuth = (Math.PI * (1.0 + Math.Sqrt(5.0)) * index) % (2.0 * Math.PI);
            var polarAngle = Math.Acos(1.0 - 2.0 * index / pointCount);
            points.Add(new SphCoord(azimuth, polarAngle, radius: null));
        }
        return points;
    }
}

/// <summary>Represents a set of spherical coordinates.
///     If Radius is given: a vector from the origin, locating a point in space.
///     If Radius is omitted or null: Direction-only coordinates (could be considered at infinite distance).
/// </summary>
public record SphCoord {
    
    public double Azimuth { get; init; }
    public double PolarAngle { get; init; }
    public double? Radius { get; init; } // Nullable.

    // TODO: Ensure Azimuth is in [0, 2*pi) and PolarAngle is in [0, pi]. >>> but HOW????? <<< 
    /// <summary>Default constructor, usage:
    ///     sc = new SphCoord {Azimuth=az, PolarAngle=pa, Radius=r};</summary>
    public SphCoord() {}

    /// <summary>Constructor from 2 or 3 (with radius) doubles.</summary>
    /// <param name="azimuth">In radians, within [0, 2*pi).</param>
    /// <param name="polarAngle">In radians, from zero at +z pole to pi at -z pole.</param>
    /// <param name="radius">Normally in meters from origin.
    ///     Set to null if omitted (object represents .</param>
    public SphCoord(double azimuth, double polarAngle, double? radius = null) {
        Azimuth = azimuth; PolarAngle = polarAngle; Radius = radius; 
    }

    /// <summary>Constructor from Point3D object (yields vector SphCoord, with valid Radius).</summary>
    /// <param name="point">A Cartesian Point3D object, here converted to spherical coordinates.</param>
    public SphCoord(Point3D point) : this(point.VectorFrom(new Point3D(0.0, 0.0, 0.0))) {
    }

    /// <summary>Constructor from Vector3D object  (yields vector SphCoord, with valid Radius).</summary>
    /// <param name="vector">A Cartesian Vector3D object, here converted to spherical coordinates.</param>
    public SphCoord(Vector3D vector) {
        Azimuth = Math.Atan2(vector.Dy, vector.Dx);
        var radius = vector.Length;
        PolarAngle = Math.Acos(vector.Dz / radius); // Cannot use property Radius here (it's nullable).
        Radius = radius;
    }
    
    

}



