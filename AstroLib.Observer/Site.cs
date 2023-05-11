// ReSharper disable InconsistentNaming
// ReSharper disable InvalidXmlDocComment

using AstroLib.Core;

namespace AstroLib.Observer;

/// <summary>Represents a specific location on or just above the earth's surface.
/// Most often specifies an observing location, i.e., primary mirror center or pivot point of the
/// telescope doing the observing. </summary>
public class Site {
    /// <summary> User's name for this site.</summary>
    public string? Name { get; set; }
    /// <summary> IAU Telescope ID for this site (e.g., "V16").</summary>
    public string? IAU { get; set; }
    
    /// <summary> Longitude (WGS84) of the observing location, in degrees.
    /// East positive, west negative.</summary>
    public double Longitude { get; init; }
    /// <summary> Longitude in hours relative to Greenwich (zero) meridian. Read-only.
    /// Respects International Date Line if OffsetFromUTC is supplied.</summary>
    public double LongitudeHours { get; private set; }
    /// <summary>Site's offset in hours from UTC, if supplied in constructor.</summary>
    public double? OffsetFromUTC { get; init; }    
    
    /// <summary> Latitude (WGS84) of the observing location, in degrees.</summary>
    public double Latitude { get; init; }
    /// <summary> Latitude (WGS84 ellipsoid) of the observing location, in meters.</summary>
    public double Elevation { get; init; }


    /// <summary>CONSTRUCTOR.</summary>
    /// <param name="longitudeDeg">Longitude (WGS84) of the observing location, in degrees.
    /// East positive, west negative.</param>
    /// <param name="latitudeDeg">Latitude (WGS84) of the observing location, in degrees.</param>
    /// <param name="elevationMeters">Latitude (WGS84 ellipsoid) of the observing location,
    /// in meters.</param>
    /// <param name="offsetFromUTC"> Offset in hours of this site from UTC time. Generally an integer,
    /// negative for the Americas, positive for Europe, Africa, Asia, ANZ. If omitted, it is derived
    /// from the site's Longitude, which will cause errors only where the International Date Line
    /// deviates from 180 Longitude (rare; e.g., Midway Island, Bering Strait, Tonga, Samoa).
    public Site(double longitudeDeg, double latitudeDeg, double elevationMeters, double? offsetFromUTC) {
        Longitude = longitudeDeg;
        Latitude = latitudeDeg;
        Elevation = elevationMeters;
        LongitudeHours = GetLongitudeInHours(offsetFromUTC);
    }

    private double GetLongitudeInHours(double? offsetFromUTC) {
        if (offsetFromUTC is null or >= 0.0) {
            return AstroMath.Wrap(Longitude, -180, 180) / 15.0;
        }
        else {
            // Site is west of the Greenwich meridian but east of the International Date Line (Americas):
            return (AstroMath.Wrap(Longitude, 0, 360) - 360.0) / 15.0;
        }
    }
}