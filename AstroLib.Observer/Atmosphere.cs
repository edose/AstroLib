// ReSharper disable InconsistentNaming

using AstroLib.Core;

namespace AstroLib.Observer; 

/// <summary>Properties of the atmosphere at the Site,
/// through which the telescope observes the sky.</summary>
// TODO: Use properties' summer and winter values, to estimate the current night's (midnight) properties.
public class Atmosphere {
    public double? ExtinctionSG { get; set; }
    public double  ExtinctionSR { get; set; }
    public double? ExtinctionSI { get; set; }
    public double? ExtinctionSZ { get; set; }
    public double RelHumidity { get; set; }
    public double Temperature { get; set; }
    
    // TODO: Add the summer/winter interpolation facilities.
    // The following are used only if class Atmosphere is being asked to interpolate between
    // typical summer and winter midnight values.
    // public double? ExtinctionSummerSG { get; set; }
    // public double? ExtinctionWinterSG { get; set; }
    // public double? ExtinctionSummerSR { get; set; }
    // public double? ExtinctionWinterSR { get; set; }
    // public double? ExtinctionSummerSI { get; set; }
    // public double? ExtinctionWinterSI { get; set; }
    // public double? ExtinctionSummerSZ { get; set; }
    // public double? ExtinctionWinterSZ { get; set; }
    // public double? RelHumiditySummer { get; set; }
    // public double? RelHumidityWinter { get; set; }
    // public double? TemperatureSummer { get; set; }
    // public double? TemperatureWinter { get; set; }
    // public int? CurrentDate      { get; set; }
    // public int? ColdestDateMonth { get; set; }
    // public int? ColdestDateDay   { get; set; }

    /// <summary> CONSTRUCTOR, from user-specified extinction, relative humidity, and temperature.</summary>
    /// <param name="extinctionSR">Zenith extinction through Sloan R filter, in magnitudes,
    /// e.g., 0.15.</param>
    /// <param name="relHumidity">Relative humidity, as fraction in range [0,1]
    /// or as percentage in range (1,100].</param>
    /// <param name="temperature">Ambient temperature, in degrees Celsius.</param>
    public Atmosphere(double extinctionSR = 0.15, double relHumidity = 0.5, double temperature = 10.0) {
        ExtinctionSR = extinctionSR;
        RelHumidity = AstroMath.Clamp((relHumidity <= 1.0) ? 
            relHumidity : relHumidity / 100.0, 0, 1.0);
        Temperature = temperature;
    }
}