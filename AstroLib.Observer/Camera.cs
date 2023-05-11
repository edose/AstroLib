// ReSharper disable InconsistentNaming
// ReSharper disable CollectionNeverQueried.Global
// ReSharper disable UnusedAutoPropertyAccessor.Global

namespace AstroLib.Observer; 

/// <summary>Represents a specific telescope, including characteristics of the optics and mount, but
/// typically not including characteristics of the filters, guiding, or camera (which are handled by
/// the Camera class). In this version, the Telescope class does not handle: pier side or flip, or
/// polar or horizon alignment of the mount.</summary>
public class Camera {
    
    /// <summary> Name of the Camera, e.g., "SBIG AC 4040".</summary>
    public string Model { get; init; }
    /// <summary> Number of pixels in sensor's X (image width) direction.</summary>
    public int XPixels { get; init; }
    /// <summary> Number of pixels in sensor's Y (image height) direction.</summary>
    public int YPixels { get; init; }
    /// <summary> Size of each pixel, in microns. Square pixels only.</summary>
    public double PixelSize { get; init; }
    /// <summary> Estimate of number of ADU counts per pixel indicating the onset of saturation.
    /// Typically, pixels having ADU count above SaturationCounts are considered invalid data.</summary>
    public double SaturationCounts { get; init; }
    /// <summary> Percent decrease in ADU count at image corner vs at image center for the same light
    /// flux. Typically in range 10 to 40. Used mostly to correct saturation counts, or for rough
    /// flat calibration of images should a flat calibration image not be available.</summary>
    public double VignettingPctAtCorner { get; init; }
    /// <summary> A list of filters (as Filter objects) available in the filter wheel.
    /// Added after constructing Camera object, via one or more calls to AddFilter()
    /// and/or AddFilters().</summary>
    public List<Filter> Filters { get; private set; }
    /// <summary> A list of transforms (as Transform objects) available.
    /// Added after constructing Camera object, via one or more calls to AddTransform()
    /// and/or AddTransforms().</summary>
    public List<Transform> Transforms { get; private set; }
    
    /// <summary> CONSTRUCTOR, normal case.</summary>
    /// <param name="model">Name of the Camera, e.g., "SBIG AC 4040".</param>
    /// <param name="xPixels">Number of pixels in sensor's X (image width) direction.</param>
    /// <param name="yPixels">Number of pixels in sensor's Y (image height) direction.</param>
    /// <param name="pixelSize">Size of each pixel, in microns. Square pixels only.</param>
    /// <param name="saturationCounts">Estimate of number of ADU counts per pixel indicating the onset of saturation.
    /// Typically, pixels having ADU count above SaturationCounts are considered invalid data.</param>
    /// <param name="vignettingPctAtCorner">Percent decrease in ADU count at image corner vs at image center for the same light
    /// flux. Typically in range 10 to 40. Used mostly to correct saturation counts, or for rough
    /// flat calibration of images should a flat calibration image not be available.</param>
    public Camera(string model, int xPixels, int yPixels, double pixelSize, double saturationCounts, 
        double vignettingPctAtCorner) {
        Model = model;
        XPixels = xPixels;
        YPixels = yPixels;
        PixelSize = pixelSize;
        SaturationCounts = saturationCounts;
        VignettingPctAtCorner = vignettingPctAtCorner;
        Filters = new List<Filter>();
        Transforms = new List<Transform>();
    }

    /// <summary> Add a Filter object to the list of filters available to the camera. </summary>
    /// <param name="filter">Filter object to add to the available list.</param>
    public void AddFilter(Filter filter) { Filters.Add(filter); }
    
    /// <summary> Add a list of Filter objects (zero, one, or more) to the list of filters
    /// available to the camera. </summary>
    /// <param name="filters">List of Filter objects to add to the available list.</param>
    public void AddFilters(IEnumerable<Filter> filters) { Filters.AddRange(filters); }
    
    /// <summary> Add a Transform object to the list of transforms suitable to the camera. </summary>
    /// <param name="transform">Transform object to add to the available list.</param>
    public void AddTransform(Transform transform) { Transforms.Add(transform); }
    
    /// <summary> Add a list of Transform objects (zero, one, or more) to the list of transforms
    /// suitable to the camera. </summary>
    /// <param name="transforms">List of Transform objects to add to the available list.</param>
    public void AddTransforms(IEnumerable<Transform> transforms) { Transforms.AddRange(transforms); }
}

// TODO: Move this to new place?: AstroLib.Core, "Things" folder?
/// <summary> Represents a physical filter (piece of glass mounted in front of the camera).
/// Not to be confused with a passband, even with exactly the same name.</summary>
public class Filter {
    public string ID { get; init; }
}