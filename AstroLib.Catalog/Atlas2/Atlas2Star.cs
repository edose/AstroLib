// ReSharper disable InconsistentNaming

using AstroLib.Core;

namespace AstroLib.Catalog.Atlas2; 

/// <summary> Represents one star from the Atlas refcat2 catalog. </summary>
public class Atlas2Star {
    
    // Identifier, e.g., 
    public string? ID { get; internal set; }
    
    // POSITION DATA at catalog epoch only. Must be corrected for proper motion before use:
    /// <summary>Right Ascension in degrees (J2000) epoch 2015.5, Gaia DR2 </summary>
    public double RaGaia { get; init; }
    /// <summary>Declination in degrees (J2000) epoch 2015.5, Gaia DR2 </summary>
    public double DecGaia { get; init; }
    
    // POSITION DATA at user's date, as corrected for proper motion:
    /// <summary>User date, for proper motion updating of RA and Declination from catalog epoch.</summary>
    public DateTime? OfDate = null;
    /// <summary>Right Ascension in degrees (J2000), at user date. Null on first construction.</summary>
    public double? RaOfDate = null;
    /// <summary>Declination in degrees (J2000), at user date. Null on first construction.</summary>
    public double? DecOfDate = null;
    
    // PARALLAX DATA (estimate of distance to star in favorable cases):
    /// <summary>Parallax in arcseconds, Gaia DR2</summary>
    public double Parallax { get; init; }
    /// <summary>Uncertainty of parallax in arcseconds, Gaia DR2</summary>
    public double uParallax { get; init; }
    
    // PROPER MOTIONS:
    /// <summary>ProperMotion in Right Ascension, in arcseconds/year, Gaia DR2</summary>
    public double PM_Ra { get; init; }
    /// <summary>ProperMotion in Declination, in arcseconds/year, Gaia DR2</summary>
    public double PM_Dec { get; init; }
    
    // NEIGHBORING FLUXES by radius around star:
    /// <summary>Radius where cumulative flux (in G) exceeds 0.1 * this star's flux, arcseconds.</summary>
    public double? RP1 { get; init; }
    /// <summary>Radius where cumulative flux (in G) exceeds 1 * this star's flux, arcseconds.</summary>
    public double? R1 { get; init; }
    /// <summary>Radius where cumulative flux (in G) exceeds 10 * this star's flux, arcseconds.</summary>
    public double? R10 { get; init; }
    
    // MAGNITUDE and COLOR DATA:
    /// <summary>Sloan g (SG) magnitude, Pan-STARRS</summary>
    public double SG { get; init; }
    /// <summary>Uncertainty in Sloan g (SG) magnitude, Pan-STARRS</summary>
    public double uSG { get; init; }
    /// <summary>Sloan r (SR) magnitude, Pan-STARRS</summary>
    public double SR { get; init; }
    /// <summary>Uncertainty in Sloan r (SR) magnitude, Pan-STARRS</summary>
    public double uSR { get; init; }
    /// <summary>Sloan i (SI) magnitude, Pan-STARRS</summary>
    public double SI { get; init; }
    /// <summary>Uncertainty in Sloan i (SI) magnitude, Pan-STARRS</summary>
    public double uSI { get; init; }
    /// <summary>Sloan z (SZ) magnitude, Pan-STARRS</summary>
    public double SZ { get; init; }
    /// <summary>Uncertainty in Sloan z (SZ) magnitude, Pan-STARRS</summary>
    public double uSZ { get; init; }
    /// <summary>r-i color, derived from SR and SI magnitudes, Pan-STARRS</summary>
    public double RIColor { get; init; }
    
    // SCREENING FLAGS:
    /// <summary>Variable flag. True -> suspicion of variability, thus unsuitable as comp star.</summary>
    public bool Variable { get; init; }
    /// <summary>NotAvailable flag. True -> some data absent, likely unsuitable as comp star.</summary>
    public bool NotAvailable { get; init; }
    /// <summary>Duplicate flag. True -> probably unsuitable as comp star.</summary>
    public bool Duplicate { get; init; }
    
    public Atlas2Star(string catalogTextLine) {
        ID = null;
        var fields = catalogTextLine.Split(',');
        RaGaia     = double.Parse(fields[0]) * 0.000_000_01;
        DecGaia    = double.Parse(fields[1]) * 0.000_000_01;
        RaOfDate   = null;
        DecOfDate  = null;
        Parallax   = double.Parse(fields[2]) * 0.000_01;  
        uParallax  = double.Parse(fields[3]) * 0.000_01;
        PM_Ra      = double.Parse(fields[4]) * 0.000_01;
        PM_Dec     = double.Parse(fields[6]) * 0.000_01;
        RP1 = (fields[18] == "999") ? null : double.Parse(fields[18]) * 0.1;
        R1  = (fields[19] == "999") ? null : double.Parse(fields[19]) * 0.1;
        R10 = (fields[20] == "999") ? null : double.Parse(fields[20]) * 0.1;
        SG  = double.Parse(fields[21]) * 0.001;
        uSG = double.Parse(fields[22]) * 0.001;
        SR  = double.Parse(fields[25]) * 0.001;
        uSR = double.Parse(fields[26]) * 0.001;
        SI  = double.Parse(fields[29]) * 0.001;
        uSI = double.Parse(fields[30]) * 0.001;
        SZ  = double.Parse(fields[33]) * 0.001;
        uSZ = double.Parse(fields[34]) * 0.001;
        RIColor = SR - SI;
        int flag     = int.Parse(fields[16]);
        Variable     = (flag % 2) == 1;
        NotAvailable = ((flag / 2) % 2) == 1;  // not of much importance in Gaia DR2/Atlas2
        Duplicate    = ((flag / 4) % 2) == 1;
    }

    /// <summary>Calculate years passed from Catalog Epoch to given date. Typically used to pre-compute
    /// a date once, for a large number of comp stars.</summary>
    /// <param name="date">Date for which RA and Dec of date are wanted.</param>
    /// <returns></returns>
    public static double YearsFromCatalogEpoch(DateTime date) {
        return (date - Atlas2Catalog.CatalogEpoch).Days / Constants.SolarDaysPerJulianYear;
    }
    
    /// <summary>Calculate and store RA and Dec for date given.
    /// This overload is rarely used; normal application uses the (double) overload with pre-computed
    /// yearsFromCatalogEpoch.</summary>
    /// <param name="date">Date for which RA and Dec of date are wanted.</param>
    public void UpdateRaDecForDate(DateTime date) {
        UpdateRaDecForDate(YearsFromCatalogEpoch(date), date);
    }

    /// <summary>Calculate and store RA and Dec for number of years passed since catalog epoch.
    /// The preferred overload, requiring the pre-computed time span.</summary>
    /// <param name="yearsFromCatalogEpoch">Years passed since catalog epoch
    /// (before epoch -> negative).</param>
    /// <param name="ofDate">User date, needed to update the OfDate property.</param>
    public void UpdateRaDecForDate(double yearsFromCatalogEpoch, DateTime ofDate) {
        OfDate = ofDate;
        var raOfDateRaw  = RaGaia  + (yearsFromCatalogEpoch * PM_Ra)  / 3600.0;
        RaOfDate  = AstroMath.EuclidianModulo(RaGaia, 360);
        DecOfDate = DecGaia + (yearsFromCatalogEpoch * PM_Dec) / 3600.0; 
    }
    
    public override string ToString() {
        return $"Atlas2Star instance at Gaia RA {RaGaia:F3}{Text.Degree}, " +
               $"Dec {DecGaia:F3}{Text.Degree}.";
    }
}

/// <summary> Represents a collection of stars from the Atlas refcat2 catalog. </summary>
public class Atlas2Stars {
    
    /// <summary>The backing property for class Atlas2Stars; a list of Atlas2Star objects.
    /// Not for user access; users may access individual star objects via the class indexer.</summary>
    public List<Atlas2Star> Stars { get; private set; }
    
    public int Count => Stars.Count;

    /// <summary>Constructor most often used, from catalog and range of RA and Dec.</summary>
    /// <param name="catalog">An Atlas2Catalog object.</param>
    /// <param name="raLow">Lowest (i.e., westernmost) Right Ascension to include.</param>
    /// <param name="raHigh">Highest (i.e., easternmost) Right Ascension to include.</param>
    /// <param name="decLow"></param>
    /// <param name="decHigh"></param>
    /// <param name="date"></param>
    public Atlas2Stars(Atlas2Catalog catalog, double raLow, double raHigh, double decLow, double decHigh, 
        DateTime? date=null) {
        var ofDate = date ?? DateTime.UtcNow;
        
        // For each of RA and Dec, make list of integers representing square-degree files to read,
        // respecting any RA zero-crossing:
        int intRaFirst, intRaLast;
        List<int> intAllRa;
        if (raLow < raHigh) {
            // RA does not cross zero:
            intRaFirst = (int) Math.Floor(raLow);
            intRaLast  = (int) Math.Floor(raHigh);
            intAllRa = Enumerable.Range(intRaFirst, intRaLast - intRaFirst + 1).ToList();
        }
        else {
            // RA does cross zero:
            intRaFirst = (int) Math.Floor(raLow);
            intRaLast  = (int) Math.Floor(raHigh);
            var intRaWest = Enumerable.Range(intRaFirst, 359 - intRaFirst + 1);
            var intRaEast = Enumerable.Range(0, intRaLast + 1);
            intAllRa = intRaWest.Concat(intRaEast).ToList();
        }
        var intDecFirst = (int) Math.Floor(decLow);
        var intDecLast  = (int) Math.Floor(decHigh);
        var intAllDec = Enumerable.Range(intDecFirst, intDecLast - intDecFirst + 1).ToList();

        // Build a master Atlas2Stars object by calling Atlas2Stars(catalog, ra, dec) for each (ra,dec):
        var masterAtlas2Stars = new Atlas2Stars(new List<Atlas2Star>());  // begin w/o stars.
        foreach(var ra in intAllRa) {
            foreach (var dec in intAllDec) {
                masterAtlas2Stars.AddStarsFrom(new Atlas2Stars(catalog, ra, dec));
            }
        }
        
        // Update date for all stars:
        masterAtlas2Stars.UpdateAllRaDecForDate(ofDate);
        
        // Screen for exact RA and Dec limits (respecting any RA zero crossing):
        List<Atlas2Star> croppedStars;
        if (raLow < raHigh) {
            // RA does not cross zero:
            croppedStars = masterAtlas2Stars.Stars
                .Where(star => star.RaOfDate  >= raLow)
                .Where(star => star.RaOfDate  <= raHigh)
                .Where(star => star.DecOfDate >= decLow)
                .Where(star => star.DecOfDate <= decHigh)
                .ToList();
        }
        else {
            // RA does cross zero:
            croppedStars = masterAtlas2Stars.Stars
                .Where(star => (star.RaOfDate  >= raLow) || (star.RaOfDate  <= raHigh))
                .Where(star => star.DecOfDate >= decLow)
                .Where(star => star.DecOfDate <= decHigh)
                .ToList();
        }
        Stars = croppedStars;
    }
    
    /// <summary>Constructor, reading and parsing one catalog text file.</summary>
    /// <param name="catalog">An Atlas2Catalog object.</param>
    /// <param name="ra">The integer (floor) Right Ascension that defines the text file to read.</param>
    /// <param name="dec">The integer (floor) Declination that defines the text file to read.</param>
    public Atlas2Stars(Atlas2Catalog catalog, int ra, int dec) {
        var stars = new List<Atlas2Star>();
        foreach (var subdirectory in catalog.CatalogSubdirectoryNamesPresent) {
            var index = 1;
            var lines = catalog.readLinesFromOneAtlasFile(subdirectory, ra, dec);
            foreach (var line in lines) {
                var star = new Atlas2Star(line);
                star.ID = $"{ra:000}{dec:+00;-00}" +
                          $"{Atlas2Catalog.SubdirectoryCodes[subdirectory]}{index:00000}";
                index++;
                stars.Add(star);
            }
        }
        Stars = stars;
    }

    /// <summary>Constructor, from plain list of Atlas2Star objects; mostly for testing.</summary>
    /// <param name="starList">A List of Atlas2Star objects.</param>
    public Atlas2Stars(List<Atlas2Star> starList) {
        Stars = starList;
    }
    
    /// <summary>Indexer, pointing to i-th star (Atlas2Star object); mostly for testing.</summary>
    /// <param name="i">Zero-based index.</param>
    public Atlas2Star this[int i] => Stars[i];

    /// <summary>Add all stars from another Atlas2Stars object to this Atlas2Stars object.</summary>
    /// <param name="otherStars">Atlas2Stars object to add to this one.</param>
    /// <returns>Total count of stars in newly expanded object.</returns>
    public int AddStarsFrom(Atlas2Stars otherStars) {
        Stars.AddRange(otherStars.Stars);
        return Stars.Count;
    }

    /// <summary>Update all stars for the given user date.</summary>
    /// <param name="date"></param>
    public void UpdateAllRaDecForDate(DateTime date) {
        double yearsFromCatalogEpoch = Atlas2Star.YearsFromCatalogEpoch(date);
        foreach (var star in Stars) {
            star.UpdateRaDecForDate(yearsFromCatalogEpoch, date);
        }
    }
    
    /// <summary>Remove all stars for which either the Variable flag or the Duplicate flag is set.</summary>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int RemoveAllFlagged() {
        Stars = Stars
            .Where(star => star.Variable     is false)
            .Where(star => star.Duplicate    is false)
            .ToList();
        return Stars.Count;
    }
    
    /// <summary>Remove each star for which catalog indicates any nearby flux.</summary>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int RemoveAllWithNearbyFlux() {
        Stars = Stars
            .Where(star => star.RP1 == null)
            .Where(star => star.R1  == null)
            .Where(star => star.R10 == null)
            .ToList();
        return Stars.Count;
    }

    /// <summary>Keep only stars in specified range of Sloan g magnitudes, and with a specified maximum
    /// listed uncertainty in Sloan g. Pass in null to any parameter to skip that screen.</summary>
    /// <param name="minMag">Minimum allowable Sloan g magnitude, or null to skip this screen.</param>
    /// <param name="maxMag">Maximum allowable Sloan g magnitude, or null to skip this screen.</param>
    /// <param name="maxUMag">Maximum allowable uncertainty in Sloan g magnitude, or null to skip.</param>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int SelectOnSG(double? minMag = null, double? maxMag = null, double? maxUMag = null) {
        Stars = Stars
            .Where(star => (minMag  == null || star.SG  >= minMag))
            .Where(star => (maxMag  == null || star.SG  <= maxMag))
            .Where(star => (maxUMag == null || star.uSG <= maxUMag))
            .ToList();
        return Stars.Count;
    }
    
    /// <summary>Keep only stars in specified range of Sloan r magnitudes, and with a specified maximum
    /// listed uncertainty in Sloan r. Pass in null to any parameter to skip that screen.</summary>
    /// <param name="minMag">Minimum allowable Sloan r magnitude, or null to skip this screen.</param>
    /// <param name="maxMag">Maximum allowable Sloan r magnitude, or null to skip this screen.</param>
    /// <param name="maxUMag">Maximum allowable uncertainty in Sloan r magnitude, or null to skip.</param>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int SelectOnSR(double? minMag = null, double? maxMag = null, double? maxUMag = null) {
        Stars = Stars
            .Where(star => (minMag  == null || star.SR  >= minMag))
            .Where(star => (maxMag  == null || star.SR  <= maxMag))
            .Where(star => (maxUMag == null || star.uSR <= maxUMag))
            .ToList();
        return Stars.Count;
    }
    
    /// <summary>Keep only stars in specified range of Sloan i magnitudes, and with a specified maximum
    /// listed uncertainty in Sloan i. Pass in null to any parameter to skip that screen.</summary>
    /// <param name="minMag">Minimum allowable Sloan i magnitude, or null to skip this screen.</param>
    /// <param name="maxMag">Maximum allowable Sloan i magnitude, or null to skip this screen.</param>
    /// <param name="maxUMag">Maximum allowable uncertainty in Sloan i magnitude, or null to skip.</param>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int SelectOnSI(double? minMag = null, double? maxMag = null, double? maxUMag = null) {
        Stars = Stars
            .Where(star => (minMag  == null || star.SI  >= minMag))
            .Where(star => (maxMag  == null || star.SI  <= maxMag))
            .Where(star => (maxUMag == null || star.uSI <= maxUMag))
            .ToList();
        return Stars.Count;
    }
    
    /// <summary>Keep only stars in specified range of Sloan z magnitudes, and with a specified maximum
    /// listed uncertainty in Sloan z. Pass in null to any parameter to skip that screen.</summary>
    /// <param name="minMag">Minimum allowable Sloan z magnitude, or null to skip this screen.</param>
    /// <param name="maxMag">Maximum allowable Sloan z magnitude, or null to skip this screen.</param>
    /// <param name="maxUMag">Maximum allowable uncertainty in Sloan z magnitude, or null to skip.</param>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int SelectOnSZ(double? minMag = null, double? maxMag = null, double? maxUMag = null) {
        Stars = Stars
            .Where(star => (minMag  == null || star.SZ  >= minMag))
            .Where(star => (maxMag  == null || star.SZ  <= maxMag))
            .Where(star => (maxUMag == null || star.uSZ <= maxUMag))
            .ToList();
        return Stars.Count;
    }
    
    /// <summary>Keep only stars in specified range of Sloan r-i color.
    /// Pass in null to either parameter to skip that screen.</summary>
    /// <param name="minRIColor">Minimum allowable Sloan r-i color, or null to skip this screen.</param>
    /// <param name="maxRIColor">Maximum allowable Sloan r-i color, or null to skip this screen.</param>
    /// <returns>Total count of stars remaining in newly screened object.</returns>
    public int SelectOnRIColor(double? minRIColor, double? maxRIColor) {
        Stars = Stars
            .Where(star => (minRIColor == null || star.RIColor >= minRIColor))
            .Where(star => (maxRIColor == null || star.RIColor <= maxRIColor))
            .ToList();
        return Stars.Count;
    }

    public override string ToString() {
        return $"Atlas2Stars: {Stars.Count} star records.";
    }
} 