namespace AstroLib.Core; 

/// <summary> Contains characters and strings commonly needed when writing Astronomy software in C#.
/// NB: We express all single characters as strings, to ease composing of user strings.
/// NB: We exclude for now: Greek letters.</summary>
// Unicode definitions retrieved 2023-04 via https://www.unicode.org/charts/.
public static class Text {

    // Format strings:
    public const string Iso8601T = "yyyy-MM-ddTHH:mm:ss.fff";  // machine-parsable; here, to milliseconds
    public const string Iso8601  = "yyyy-MM-dd HH:mm:ss.fff";  // human-preferred; here, to milliseconds

    // Modifiers and units:
    public const string Degree = "\u00B0";
    public const string Micro = "\u00B5";  // ~ Greek small micron.
    public const string Angstrom = "\u212b";
    
    // Operators & relations:
    public const string PlusMinus = "\u00B1";
    public const string MultiplyBy = "\u00D7";
    public const string DivideBy = "\u00F7";
    public const string Infinity = "\u221E";
    public const string AlmostEqualTo = "\u2248";
    public const string ApproxEqualTo = "\u2245";
    public const string NotEqualTo = "\u2260";
    public const string IdenticalTo = "\u2261";
    public const string MuchLessThan = "\u00AB";
    public const string MuchGreaterThan = "\u00BB";    
    public const string LessThanOrEqualTo = "\u2264";
    public const string GreaterThanOrEqualTo = "\u2265";
    public const string Therefore = "\u2234";
    public const string Because = "\u2235";

    // Astronomy glyphs:
    public const string Telescope = "\u0001F52D";
    public const string SatelliteDish = "\u0001F4E1";
    
    public const string WhiteStar = "\u2606";
    public const string BlackStar = "\u2605";
    
    public const string Comet = "\u2604";
    public const string RingedPlanet = "\u0001FA90";
    
    public const string WhiteSunWithRays = "\u263C";
    public const string BlackSunWithRays = "\u2600";
    public const string SunWithFace = "\u0001F31E";    
    
    public const string MoonCrescent = "\U0001F319";
    public const string FirstQuarterMoonRound = "\u0001F313";
    public const string FirstQuarterMoonCrescent = "\u263D";
    public const string LastQuarterMoonRound = "\u0001F317";
    public const string LastQuarterMoonCrescent = "\u263E";

    // Weather:
    public const string Cloud = "\u2601";
    public const string CloudWithRain = "\u0001F327";
    public const string Umbrella = "\u2602";
    public const string UmbrellaWithRain = "\u2614";
    public const string SunBehindCloud = "\u26C5";
    public const string SunWithSmallCloud = "\u0001F324";
    
    // Other glyphs:
    public const string CancellationX = "\u0001F5D9";
    public const string CrossMark = "\u274C";
    public const string CheckMark = "\u2713";
    public const string OpenBook = "\u0001F56E";
    public const string OpenBookWithLines = "\u0001F4D6";
    public const string Satellite = "\u0001F6F0";

    // For use in text blocks:
    public const string BulletPoint = "\u2022";
    public const string Ellipsis = "\u2026";
    public const string NoBreakSpace = "\u00A0";  // = NBSP
    public const string EmDash = "\u2014";
    public const string Copyright = "\u00A9";
    public const string Trademark = "\u2122";
}

