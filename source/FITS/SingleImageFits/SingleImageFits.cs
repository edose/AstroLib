using System;
using System.Diagnostics.CodeAnalysis;

namespace AstroLib.Fits.SingleImageFits;

/// <summary>Represents a Single Image FITS file, as defined in FITS standard.
/// One HDU only (no extensions); 2-dimensional (monochrome image) data only.</summary>
public class SingleImageFits
{ 
    public HduHeader Header { get; private set; }
    public HduData Data { get; private set; }    
    public FitsDataType DataType { get; private set; }
    
    /// <summary>Constructor, from full path of FITS file.</summary>
    /// <param name="fullPath">Full path to FITS file to be read, e.g. 'C:/Astro/FITS/fits01.fts'.</param>
    public SingleImageFits(string fullPath)
    {
        (Header, Data) = ReadAndParseSif(fullPath);
        
        // TODO: make & expose (lots of) convenience properties.
    }

    /// <summary>Read and parse Single Image FITS file, keeping file open as briefly as possible.</summary>
    /// <param name="fullPath">Full path to FITS file to be read.</param>
    /// <returns>2-Tuple of a parsed HduHeader object and a parsed HduData object.</returns>
    private static (HduHeader hduHeader, HduData? hduData) ReadAndParseSif(string fullPath)
    {
        using var fs = new FileStream(fullPath, FileMode.Open);
        fs.Position = 0;
        
        var hduHeader = new HduHeader(fs);
        bool canReadData = hduHeader.MandatoryKeywordsPresent &&
                           (hduHeader.ValueRead["NAXIS"].ValueInteger == 2) &&
                           (hduHeader.ValueRead["NAXIS1"].ValueInteger > 0) &&
                           (hduHeader.ValueRead["NAXIS2"].ValueInteger > 0);
        var hduData = (canReadData) ? new HduData(fs, hduHeader) : null;
        return (hduHeader, hduData);
    }
}