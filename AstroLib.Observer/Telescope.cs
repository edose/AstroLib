namespace AstroLib.Observer; 

/// <summary>Represents a specific telescope, including characteristics of the optics and mount, but
/// typically not including characteristics of the filters, guiding, or camera (which are handled by
/// the Camera class). In this version, the Telescope class does not handle: pier side or flip, or
/// polar or horizon alignment of the mount.</summary>
public class Telescope {
    /// <summary> Name of the telescope mount, e.g., "PlaneWave L-500".</summary>
    public string MountModel { get; init; }
    /// <summary> Name of the Optical Tube Assembly (telescope optics proper),
    /// e.g., "Celestron C14".</summary>
    public string OtaModel { get; init; }
    
    /// <summary> Aperture of the OTA, in meters,
    /// e.g., 0.350 for OTA having a 14-inch primary mirror.</summary>
    public double Aperture { get; init; }
    /// <summary> Focal length of the OTA including any focal reducer, in meters,
    /// e.g. 2.800 for OTA having f/8 mirror and 0.35 meter aperture.</summary>
    public double FocalLength { get; init; }
    
    /// <summary>CONSTRUCTOR for normal case.</summary>
    /// <param name="mountModel">Name of the telescope mount, e.g., "PlaneWave L-500".</param>
    /// <param name="otaModel">Name of the Optical Tube Assembly (telescope optics proper),
    /// e.g., "Celestron C14".</param>
    /// <param name="focalLength">Focal length of the OTA, in meters,
    /// e.g., 350 for OTA having a 14-inch primary mirror.</param>
    public Telescope(string mountModel, string otaModel, double focalLength) {
        MountModel = mountModel;
        OtaModel = otaModel;
        FocalLength = focalLength;
    }
}