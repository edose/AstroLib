// ReSharper disable InconsistentNaming
// ReSharper disable MemberCanBePrivate.Global

namespace AstroLib.Core; 

/// <summary> Contains numerical constants commonly needed when writing Astronomy software in C#.</summary>
public static class Constants {
    
    // Dimension conversions:
    /// <summary>Nominal value, varies with earth rotation (i.e., with UT1). </summary>
    public const double SecondsPerSiderealDay = 86164.0905;
    /// <summary> By definition.</summary>
    public const double SecondsPerSolarDay = 24.0 * 3600.0;
    /// <summary> By definition.</summary>
    public const double SolarDaysPerJulianYear = 365.25;
    public const double SecondsPerJulianYear = SecondsPerSolarDay * SolarDaysPerJulianYear;
    /// <summary>Approximately 57.295780.</summary>
    public const double DegreesPerRadian = 180.0 / Math.PI;
    public const double ArcminutesPerRadian = 60.0 * DegreesPerRadian;
    /// <summary>Approximately 204265.</summary>
    public const double ArcsecondsPerRadian = 3600.0 * DegreesPerRadian;
    /// <summary>Commonly adopted standard, especially ISO and ICAO.</summary>
    public const double PascalsPerStandardAtmosphere = 101325;
    /// <summary>IAU definition of 2012.</summary>
    public const double MetersPerAU = 149_597_870_700;
    public const double MetersPerLightYear = SecondsPerJulianYear * LightSpeed;
    /// <summary>IAU definition of 2015.</summary>
    public const double MetersPerParsec = (180 * 60 * 60 / Math.PI) * MetersPerAU;

    // Astronomical and physical:
    /// <summary> By definition, in meters/second (SI).</summary>
    public const double LightSpeed = 299_792_458;
    /// <summary>Common alias for LightSpeed.</summary>
    public const double C = LightSpeed;
    /// <summary>CODATA definition of 2018, in Newton meters^2 / second^2 (SI).</summary>
    public const double GravitationalConstant = 6.67430E-11;
    /// <summary>Common alias for GravitationalConstant.</summary>
    public const double G = GravitationalConstant;
    /// <summary>Standard value for earth's gravitational constant, in m^3 / s^2 (SI).</summary>
    public const double GravitationalConstantEarth = 3.986004418E14;
    /// <summary>Alias for GravitationalConstantEarth.</summary>
    public const double GMEarth = GravitationalConstantEarth;
    /// <summary>Standard value for sun's gravitational constant, in m^3 / s^2 (SI).</summary>
    public const double GravitationalConstantSun = 1.32712440018E20;
    /// <summary>Alias for GravitationalConstantSun.</summary>
    public const double GMSun = GravitationalConstantSun; // alias
    /// <summary>Standard value of 1901 for earth gravitational acceleration at earth surface,
    /// in meter / second^2 (SI).</summary>
    public const double EarthGravity = 9.80665;
    /// <summary>Standard value for earth's mass, in kg (SI).</summary>
    public const double EarthMass = 5.972_167_87E+24;
    /// <summary>Standard value for sun's mass (one solar mass), in kg (SI).</summary>
    public const double SolarMass = 1.988_409_87E+30; 
    
    // Astrometric:
    /// <summary>Conversion from full-width-at-half-maximum (FWHM) to Gaussian standard deviation (sigma)
    /// for a circular, 2-dimensional Gaussian profile. Typically used in measuring width of a star profile
    /// in an astronomical image.</summary>
    public static readonly double FwhmPerSigma = 2.0 * Math.Sqrt(2.0 * Math.Log(2.0));
}