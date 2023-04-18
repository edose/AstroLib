// ReSharper disable InconsistentNaming
// ReSharper disable MemberCanBePrivate.Global

namespace AstroLib.Core; 

/// <summary> Contains numerical constants commonly needed when writing Astronomy software in C#.</summary>
public static class Constants {
    
    // Dimension conversions:
    public const double SecondsPerSiderealDay = 86164.0905;   // nominal value; varies w/ earth rot. (UT1).
    public const double SecondsPerSolarDay = 24.0 * 3600.0; // by definition    
    public const double SolarDaysPerJulianYear = 365.25; // by definition
    public const double SecondsPerJulianYear = SecondsPerSolarDay * SolarDaysPerJulianYear;
    public const double DegreesPerRadian = 180.0 / Math.PI;
    public const double ArcminutesPerRadian = 60.0 * DegreesPerRadian;
    public const double ArcsecondsPerRadian = 3600.0 * DegreesPerRadian;  // ca. 204265.
    public const double PascalsPerStandardAtmosphere = 101325;
    public const double MetersPerAU = 149_597_870_700; // IAU 2012 definition
    public const double MetersPerLightYear = SecondsPerJulianYear * LightSpeed;
    public const double MetersPerParsec = (180 * 60 * 60 / Math.PI) * MetersPerAU; // IAU 2015 definition

    // Astronomical and physical:
    public const double LightSpeed = 299_792_458; // meters/second (SI), by definition
    public const double C = LightSpeed;           // alias
    public const double GravitationalConstant = 6.67430E-11; // N m^2 / s^2 (SI); CODATA 2018 definition.
    public const double G = GravitationalConstant; // alias
    public const double GravitationalConstantEarth = 3.986004418E14; // m^3 / s^2 (SI); 
    public const double GMEarth = GravitationalConstantEarth; // alias
    public const double GravitationalConstantSun = 1.32712440018E20; // m^3 / s^2 (SI);
    public const double GMSun = GravitationalConstantSun; // alias
    public const double EarthGravity = 9.80665; // m / s^2 (SI), 1901 standard value.
    public const double EarthMass = 5.972_167_87E+24; // kg (SI)
    public const double SolarMass = 1.988_409_87E+30; // kg (SI)    
    
    // Astrometric:
    public static readonly double FwhmPerSigma = 2.0 * Math.Sqrt(2.0 * Math.Log(2.0));
}