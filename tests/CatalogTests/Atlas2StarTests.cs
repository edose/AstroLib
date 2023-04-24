using AstroLib.Catalog.Atlas2;
using AstroLib.Core;
using NUnit.Framework;
// ReSharper disable RedundantBoolCompare
// ReSharper disable InconsistentNaming

namespace AstroLibTests.CatalogTests; 

[TestFixture]
public class Atlas2StarTests {
    private string topPath;
    private Atlas2Catalog cat;
    private List<string> catLines;
    private const string FromTonryArticle = "28000001672,-1967818581,98,10,114,16,-1460,15,15884,1," + 
                                            "16472,10,15137,1,4729,895,2,1234,50,50,155,16657,10," + 
                                            "23,1f,15915,12,41,3f,15578,10,49,0f,15346,12,0,06," + 
                                            "0,14105,36,14105,53,13667,44";
    private string line44;
    
    [SetUp]
    public void Setup() {
        topPath = @"D:\Astro\Catalogs\ATLAS-refcat2";
        cat = new Atlas2Catalog(topPath);
        catLines = cat.readLinesFromOneAtlasFile("mag-0-16", 149, 71);
        line44 = catLines[44];
    }

    [Test]
    public void Class_Constructor_Tests() {
        // Test using Table 4 from Atlas refcat2 article (Tonry 2018):
        var fields = FromTonryArticle.Split(',');
        Assert.That(fields.Length, Is.EqualTo(44));
        var star = new Atlas2Star(FromTonryArticle);
        Assert.IsInstanceOf(typeof(Atlas2Star), star);
        Assert.That(star.ID, Is.Null);
        Assert.That(star.RaGaia,    Is.EqualTo(280.000_016_72).Within(1E-8));
        Assert.That(star.DecGaia,   Is.EqualTo(-19.678_185_81).Within(1E-8));
        Assert.That(star.OfDate,    Is.Null);
        Assert.That(star.RaOfDate,  Is.Null);
        Assert.That(star.DecOfDate, Is.Null);
        Assert.That(star.Parallax,  Is.EqualTo(0.000_98).Within(1E-8));
        Assert.That(star.uParallax, Is.EqualTo(0.000_10).Within(1E-8));
        Assert.That(star.PM_Ra,     Is.EqualTo(0.001_14).Within(1E-8));
        Assert.That(star.PM_Dec,    Is.EqualTo(-0.014_60).Within(1E-8));
        Assert.That(star.RP1,  Is.EqualTo(5.0));
        Assert.That(star.R1,   Is.EqualTo(5.0));
        Assert.That(star.R10,  Is.EqualTo(15.5));
        Assert.That(star.SG,  Is.EqualTo(16.657).Within(1E-8));
        Assert.That(star.uSG, Is.EqualTo( 0.010).Within(1E-8));
        Assert.That(star.SR,  Is.EqualTo(15.915).Within(1E-8));
        Assert.That(star.uSR, Is.EqualTo( 0.012).Within(1E-8));
        Assert.That(star.SI,  Is.EqualTo(15.578).Within(1E-8));
        Assert.That(star.uSI, Is.EqualTo( 0.010).Within(1E-8));
        Assert.That(star.SZ,  Is.EqualTo(15.346).Within(1E-8));
        Assert.That(star.uSZ, Is.EqualTo( 0.012).Within(1E-8));
        Assert.That(star.RIColor, Is.EqualTo(star.SR - star.SI).Within(1E-8));
        Assert.That(star.Variable, Is.False);
        Assert.That(star.NotAvailable, Is.True);
        Assert.That(star.Duplicate, Is.False);

        // Test using 45th line (string line[44]) from an actual Atlas2 file:
        var star45 = new Atlas2Star(line44);
        Assert.That(star45.ID, Is.Null);
        Assert.That(star45.RaGaia,    Is.EqualTo(149.221_786_41).Within(1E-8));
        Assert.That(star45.DecGaia,   Is.EqualTo( 71.084_643_86).Within(1E-8));
        Assert.That(star45.OfDate,    Is.Null);
        Assert.That(star45.RaOfDate,  Is.Null);
        Assert.That(star45.DecOfDate, Is.Null);
        Assert.That(star45.Parallax,  Is.EqualTo(0.000_76).Within(1E-8));
        Assert.That(star45.uParallax, Is.EqualTo(0.000_02).Within(1E-8));
        Assert.That(star45.PM_Ra,     Is.EqualTo(0.003_17).Within(1E-8));
        Assert.That(star45.PM_Dec,    Is.EqualTo(-0.010_30).Within(1E-8));
        Assert.That(star45.RP1,  Is.EqualTo(34.6));
        Assert.That(star45.R1,   Is.Null);
        Assert.That(star45.R10,  Is.Null);
        Assert.That(star45.SG,  Is.EqualTo(15.280).Within(1E-8));
        Assert.That(star45.uSG, Is.EqualTo( 0.008).Within(1E-8));
        Assert.That(star45.SR,  Is.EqualTo(14.911).Within(1E-8));
        Assert.That(star45.uSR, Is.EqualTo( 0.008).Within(1E-8));
        Assert.That(star45.SI,  Is.EqualTo(14.791).Within(1E-8));
        Assert.That(star45.uSI, Is.EqualTo( 0.008).Within(1E-8));
        Assert.That(star45.SZ,  Is.EqualTo(14.768).Within(1E-8));
        Assert.That(star45.uSZ, Is.EqualTo( 0.008).Within(1E-8));
        Assert.That(star45.RIColor, Is.EqualTo(star45.SR - star45.SI).Within(1E-8));
        Assert.That(star45.Variable, Is.False);
        Assert.That(star45.NotAvailable, Is.True);
        Assert.That(star45.Duplicate, Is.False);
    }

    [Test]
    public void Class_Property_Tests() {
        // Properties are tested in Constructor_Tests():
    }

    [Test]
    public void Method_YearsFromCatalogEpoch_Test() {
        // var star = new Atlas2Star(FromTonryArticle);
        var userDate = new DateTime(2022, 4, 1, 0, 0, 0);
        var expectedYearsFromCatalogEpoch = +6.75;
        Assert.That(Atlas2Star.YearsFromCatalogEpoch(userDate), 
            Is.EqualTo(expectedYearsFromCatalogEpoch).Within(0.02));
    }

    [Test]
    public void Method_UpdateRaDecForDate_Test() {
        var star = new Atlas2Star(FromTonryArticle);
        var userDate = new DateTime(2022, 4, 1, 0, 0, 0);
        var expectedYearsFromCatalogEpoch = +6.75;
        
        //.UpdateRaDecForDate(date):
        Assert.That(star.OfDate,    Is.Null);
        Assert.That(star.RaOfDate,  Is.Null);
        Assert.That(star.DecOfDate, Is.Null);
        star.UpdateRaDecForDate(userDate);
        Assert.That(star.OfDate, Is.EqualTo(userDate));
        Assert.That(star.RaOfDate,   
            Is.EqualTo(star.RaGaia  + expectedYearsFromCatalogEpoch * star.PM_Ra  / 3600.0)
                .Within(0.000_01));
        Assert.That(star.DecOfDate,  
            Is.EqualTo(star.DecGaia + expectedYearsFromCatalogEpoch * star.PM_Dec / 3600.0)
                .Within(0.000_01));
    
        //.UpdateRaDecForDate(double):
        var star2 = new Atlas2Star(FromTonryArticle);
        Assert.That(star2.OfDate,    Is.Null);
        Assert.That(star2.RaOfDate,  Is.Null);
        Assert.That(star2.DecOfDate, Is.Null);
        star2.UpdateRaDecForDate(expectedYearsFromCatalogEpoch, userDate);
        Assert.That(star.OfDate, Is.EqualTo(userDate));
        Assert.That(star2.RaOfDate,   
            Is.EqualTo(star2.RaGaia  + expectedYearsFromCatalogEpoch * star2.PM_Ra  / 3600.0)
                .Within(0.000_01));
        Assert.That(star2.DecOfDate,  
            Is.EqualTo(star2.DecGaia + expectedYearsFromCatalogEpoch * star2.PM_Dec / 3600.0)
                .Within(0.000_01));
    }

    [Test]
    public void Method_ToString_Test() {
        var star = new Atlas2Star(FromTonryArticle);
        Assert.That(star.ToString(), 
            Is.EqualTo($"Atlas2Star instance at Gaia RA 280.000{Text.Degree}, Dec -19.678{Text.Degree}."));
        Assert.That(star.ToString(), Is.EqualTo(star.ToString()));
    }
} // class Atlas2StarTests.

[TestFixture]
public class Atlas2StarsTests {

    private string topPath;
    private Atlas2Catalog cat;
    private readonly DateTime testDate = new DateTime(2023, 04, 23, 0, 0, 0);

    [SetUp]
    public void Setup() {
        topPath = @"D:\Astro\Catalogs\ATLAS-refcat2";
        cat = new Atlas2Catalog(topPath);
    }

    [Test]
    public void Class_Constructor_Tests() {
        // Constructor from catalog, one RA, one Dec:
        var stars = new Atlas2Stars(cat, 149, 71);
        Assert.IsInstanceOf(typeof(Atlas2Stars), stars);
        Assert.That(stars.Count, Is.EqualTo(872));

        // Constructor from pre-existing list of Atlas2Star objects:
        var fourLines = cat.readLinesFromOneAtlasFile("mag-0-16", 149, 71)
            .ToArray()[..4];
        var fourStarObjects = fourLines.Select(line => new Atlas2Star(line)).ToList();
        var fourStars = new Atlas2Stars(fourStarObjects);
        Assert.IsInstanceOf(typeof(Atlas2Stars), fourStars);
        Assert.That(fourStars.Count, Is.EqualTo(4));
        Assert.That(fourStars[0].OfDate, Is.Null);
        Assert.That(fourStars[3].OfDate, Is.Null);
        Assert.That(fourStars[0].RaGaia, Is.EqualTo(149.013_695_80).Within(1E-8));
        Assert.That(fourStars[3].DecGaia, Is.EqualTo(71.639_601_45).Within(1E-8));

        // Constructor from empty pre-existing list of Atlas2Star objects:
        var noStars = new Atlas2Stars(new List<Atlas2Star>());
        Assert.IsInstanceOf(typeof(Atlas2Stars), noStars);
        Assert.That(noStars.Count, Is.EqualTo(0));

        // Main constructor, from catalog, ranges in RA and Dec, and optional (recommended) user date:
        // Across RA and Dec degree borders (multiple files), RA does not cross zero:
        var starsA = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        Assert.IsInstanceOf(typeof(Atlas2Stars), starsA);
        Assert.That(starsA.Count, Is.EqualTo(639));
        CollectionAssert.AllItemsAreUnique(starsA.Stars);
        Assert.That(starsA.Stars.All(star => star.RaOfDate >= 23.75));
        Assert.That(starsA.Stars.All(star => star.RaOfDate <= 24.1));
        Assert.That(starsA.Stars.All(star => star.DecOfDate >= 33.7));
        Assert.That(starsA.Stars.All(star => star.DecOfDate <= 34.14));
        Assert.That(starsA.Stars.Any(star => star.RaOfDate < 23.8));
        Assert.That(starsA.Stars.Any(star => star.RaOfDate > 24.05));
        Assert.That(starsA.Stars.Any(star => star.DecOfDate < 33.8));
        Assert.That(starsA.Stars.Any(star => star.DecOfDate > 34.04));

        // Entirely within one RA and Dec degree file:
        var starsB = new Atlas2Stars(cat, 13.5, 13.7, -12.7, -12.4, testDate);
        Assert.IsInstanceOf(typeof(Atlas2Stars), starsB);
        Assert.That(starsB.Count, Is.EqualTo(88));
        CollectionAssert.AllItemsAreUnique(starsB.Stars);
        Assert.That(starsB.Stars.All(star => star.RaOfDate >= 13.5));
        Assert.That(starsB.Stars.All(star => star.RaOfDate <= 13.7));
        Assert.That(starsB.Stars.All(star => star.DecOfDate >= -12.7));
        Assert.That(starsB.Stars.All(star => star.DecOfDate <= -12.4));

        // Across RA and Dec degree boundaries, limits directly on boundaries, and RA crosses zero:
        var starsC = new Atlas2Stars(cat, 358, 2, 40, 43, testDate);
        Assert.IsInstanceOf(typeof(Atlas2Stars), starsC);
        Assert.That(starsC.Count, Is.EqualTo(78336));
        Assert.That(starsC.Stars.All(star => star.RaOfDate is (>= 358 and < 360) or (>= 0 and <= 2)));
        Assert.That(starsC.Stars.All(star => star.DecOfDate is (>= 40 and <= 43)));
    }

    [Test]
    public void Class_Property_Tests() {
        // Properties Count and indexer were tested in Constructor_Tests() above.
        
        // .NullCount:
        {
            var stars = new Atlas2Stars(cat, 149, 71);
            Assert.That(stars.NullCount, Is.EqualTo(stars.Count));
            var dateA = new DateTime(2023, 4, 23, 0, 0, 0);
            foreach (var star in stars.Stars) star.OfDate = dateA;
            Assert.That(stars.NullCount, Is.Zero);
            stars[2].OfDate = null;
            Assert.That(stars.NullCount, Is.EqualTo(1));
        }

        // .AllDatesValid:
        {
            // Case: no stars in Atlas2Stars object:
            var stars = new Atlas2Stars(new List<Atlas2Star>()); 
            Assert.That(stars.AllDatesValid == false);
            
            // Case: absent (all null) OfDate values: 
            stars = new Atlas2Stars(cat, 149, 71);
            Assert.That(stars.Count, Is.EqualTo(872));
            Assert.That(stars.AllDatesValid == false);
            
            // Case: a non-uniform OfDate value:
            var dateA = new DateTime(2023, 4, 23, 0, 0, 0);
            var dateB = new DateTime(2023, 4, 11, 11, 11, 11);
            foreach (var star in stars.Stars) star.OfDate = dateA;
            Assert.That(stars.AllDatesValid == true);
            stars[1].OfDate = dateB;
            Assert.That(stars.AllDatesValid == false);
            
            // Case: a null OfDate value: 
            foreach (var star in stars.Stars) star.OfDate = dateA;
            Assert.That(stars.AllDatesValid == true);
            stars[2].OfDate = null;
            Assert.That(stars.AllDatesValid == false);
        }
    }
    
    [Test]
    public void Method_AddStarsFrom_Test() {
        var starsA = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var starsACount = starsA.Count;
        var starsB = new Atlas2Stars(cat, 24.1, 24.3, 33.7, 34.14, testDate);
        var starsBCount = starsB.Count;
        Assert.That(starsACount, Is.EqualTo(639));
        Assert.That(starsBCount, Is.EqualTo(330));
        var allStarCount = starsA.AddStarsFrom(starsB);
        Assert.That(allStarCount, Is.EqualTo(starsACount + starsBCount));
        Assert.That(starsA.Stars.All(star => star.RaOfDate is (>= 23.75 and <= 24.3)));
        Assert.That(starsA.Stars.All(star => star.DecOfDate is (>= 33.7 and <= 34.14)));
    }
    
    [Test]
    public void Method_RemoveAllWithNearbyFlux_Test() {
        var starsC = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var beforeCount = starsC.Count;
        var beforeWithFluxCount = starsC.Stars
            .Where(star => (star.RP1 != null || star.R1 != null || star.R10 != null))
            .ToList().Count;
        Assert.That(beforeWithFluxCount, Is.Positive);
        var afterCount = starsC.RemoveAllWithNearbyFlux();
        var afterWithFluxCount = starsC.Stars
            .Where(star => (star.RP1 != null || star.R1 != null || star.R10 != null))
            .ToList().Count;
        Assert.That(afterWithFluxCount, Is.Zero);
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeWithFluxCount));
    }

    [Test]
    public void Method_RemoveAllFlaggedStars_Test() {
        var starsA = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var beforeCount = starsA.Count;
        var beforeFlaggedCount = starsA.Stars
            .Where(star => (star.Variable is true || star.Duplicate is true))
            .ToList().Count;
        Assert.That(beforeFlaggedCount, Is.Positive);
        var afterCount = starsA.RemoveAllFlagged();
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeFlaggedCount));
        var afterFlaggedCount = starsA.Stars
            .Where(star => (star.Variable is true || star.Duplicate is true))
            .ToList().Count;
        Assert.That(afterFlaggedCount, Is.Zero);
    }

    [Test]
    public void Method_SelectOnSG_Mag_Test() {
        // Test cropping of SG to range 13-14.5:
        {
            var starsD = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = starsD.Count;
            var beforeOutsideSGRange = starsD.Stars.Where(star => star.SG is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(beforeOutsideSGRange, Is.Positive);
            var afterCount = starsD.SelectOnSG(13, 14.5);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSGRange));
            var afterOutsideSGRange = starsD.Stars.Where(star => star.SG is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(afterOutsideSGRange, Is.Zero);
        }
        // Test cropping of SG with min of 13.5 only:
        {
            var starsD = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = starsD.Count;
            var beforeOutsideSGRange = starsD.Stars.Where(star => star.SG < 13.5)
                .ToList().Count;
            Assert.That(beforeOutsideSGRange, Is.Positive);
            var afterCount = starsD.SelectOnSG(minMag: 13.5, maxMag: null);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSGRange));
            var afterOutsideSGRange = starsD.Stars.Where((star => star.SG < 13.5))
                .ToList().Count;
            Assert.That(afterOutsideSGRange, Is.Zero);
        }
        // Test cropping of SG with max of 14.0 only:
        {
            var starsD = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = starsD.Count;
            var beforeOutsideSGRange = starsD.Stars.Where(star => star.SG > 14)
                .ToList().Count;
            Assert.That(beforeOutsideSGRange, Is.Positive);
            var afterCount = starsD.SelectOnSG(minMag: null, maxMag: 14);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSGRange));
            var afterOutsideSGRange = starsD.Stars.Where(star => star.SG > 14)
                .ToList().Count;
            Assert.That(afterOutsideSGRange, Is.Zero);
        }
    }

    [Test]
    public void Method_SelectOnSG_Uncertainty_Test() {
        var starsD = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var testUMag = 0.0105;  // 10.5 millimagnitudes.
        
        // Test cropping of SG uncertainty to max limit:
        var beforeCount = starsD.Count;
        var beforeUMagTooHigh = starsD.Stars.Where(star => star.uSG > testUMag)
            .ToList().Count;
        Assert.That(beforeUMagTooHigh, Is.Positive);
        var afterCount = starsD.SelectOnSG(maxUMag: testUMag);
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeUMagTooHigh));
        var afterUMagTooHigh = starsD.Stars.Where(star => star.uSG > testUMag)
            .ToList().Count;
        Assert.That(afterUMagTooHigh, Is.Zero);
    }
    
    [Test]
    public void Method_SelectOnSR_Mag_Test() {
        // Test cropping of SR to range 13-14.5:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSRRange = stars.Stars.Where(star => star.SR is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(beforeOutsideSRRange, Is.Positive);
            var afterCount = stars.SelectOnSR(13, 14.5);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSRRange));
            var afterOutsideSRRange = stars.Stars.Where(star => star.SR is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(afterOutsideSRRange, Is.Zero);
        }
        // Test cropping of SR with min of 13.5 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSRRange = stars.Stars.Where(star => star.SG < 13.5)
                .ToList().Count;
            Assert.That(beforeOutsideSRRange, Is.Positive);
            var afterCount = stars.SelectOnSG(minMag: 13.5, maxMag: null);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSRRange));
            var afterOutsideSRRange = stars.Stars.Where((star => star.SG < 13.5))
                .ToList().Count;
            Assert.That(afterOutsideSRRange, Is.Zero);
        }
        // Test cropping of SR with max of 14.0 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSRRange = stars.Stars.Where(star => star.SG > 14)
                .ToList().Count;
            Assert.That(beforeOutsideSRRange, Is.Positive);
            var afterCount = stars.SelectOnSG(minMag: null, maxMag: 14);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSRRange));
            var afterOutsideSRRange = stars.Stars.Where(star => star.SG > 14)
                .ToList().Count;
            Assert.That(afterOutsideSRRange, Is.Zero);
        }
    }
    
    [Test]
    public void Method_SelectOnSR_Uncertainty_Test() {
        var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var testUMag = 0.0105;  // 10.5 millimagnitudes.
        
        // Test cropping of SR uncertainty to max limit:
        var beforeCount = stars.Count;
        var beforeUMagTooHigh = stars.Stars.Where(star => star.uSR > testUMag)
            .ToList().Count;
        Assert.That(beforeUMagTooHigh, Is.Positive);
        var afterCount = stars.SelectOnSR(maxUMag: testUMag);
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeUMagTooHigh));
        var afterUMagTooHigh = stars.Stars.Where(star => star.uSR > testUMag)
            .ToList().Count;
        Assert.That(afterUMagTooHigh, Is.Zero);
    }
    
    [Test]
    public void Method_SelectOnSI_Mag_Test() {
        // var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        
        // Test cropping of SI to range 13-14.5:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSIRange = stars.Stars.Where(star => star.SI is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(beforeOutsideSIRange, Is.Positive);
            var afterCount = stars.SelectOnSI(13, 14.5);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSIRange));
            var afterOutsideSIRange = stars.Stars.Where(star => star.SI is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(afterOutsideSIRange, Is.Zero);
        }
        // Test cropping of SI with min of 13.5 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSIRange = stars.Stars.Where(star => star.SI < 13.5)
                .ToList().Count;
            Assert.That(beforeOutsideSIRange, Is.Positive);
            var afterCount = stars.SelectOnSI(minMag: 13.5, maxMag: null);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSIRange));
            var afterOutsideSIRange = stars.Stars.Where((star => star.SI < 13.5))
                .ToList().Count;
            Assert.That(afterOutsideSIRange, Is.Zero);
        }
        // Test cropping of SI with max of 14.0 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSIRange = stars.Stars.Where(star => star.SI > 14)
                .ToList().Count;
            Assert.That(beforeOutsideSIRange, Is.Positive);
            var afterCount = stars.SelectOnSI(minMag: null, maxMag: 14);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSIRange));
            var afterOutsideSIRange = stars.Stars.Where(star => star.SI > 14)
                .ToList().Count;
            Assert.That(afterOutsideSIRange, Is.Zero);
        }
    }
    
    [Test]
    public void Method_SelectOnSI_Uncertainty_Test() {
        var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var testUMag = 0.0105;  // 10.5 millimagnitudes.
        
        // Test cropping of SI uncertainty to max limit:
        var beforeCount = stars.Count;
        var beforeUMagTooHigh = stars.Stars.Where(star => star.uSI > testUMag)
            .ToList().Count;
        Assert.That(beforeUMagTooHigh, Is.Positive);
        var afterCount = stars.SelectOnSI(maxUMag: testUMag);
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeUMagTooHigh));
        var afterUMagTooHigh = stars.Stars.Where(star => star.uSI > testUMag)
            .ToList().Count;
        Assert.That(afterUMagTooHigh, Is.Zero);
    }

    [Test]
    public void Method_SelectOnSZ_Mag_Test() {
        // Test cropping of SZ to range 13-14.5:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSZRange = stars.Stars.Where(star => star.SZ is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(beforeOutsideSZRange, Is.Positive);
            var afterCount = stars.SelectOnSZ(13, 14.5);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSZRange));
            var afterOutsideSZRange = stars.Stars.Where(star => star.SZ is (< 13 or > 14.5))
                .ToList().Count;
            Assert.That(afterOutsideSZRange, Is.Zero);
        }
        // Test cropping of SZ with min of 13.5 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSZRange = stars.Stars.Where(star => star.SZ < 13.5)
                .ToList().Count;
            Assert.That(beforeOutsideSZRange, Is.Positive);
            var afterCount = stars.SelectOnSZ(minMag: 13.5, maxMag: null);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSZRange));
            var afterOutsideSZRange = stars.Stars.Where((star => star.SZ < 13.5))
                .ToList().Count;
            Assert.That(afterOutsideSZRange, Is.Zero);
        }
        // Test cropping of SZ with max of 14.0 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideSZRange = stars.Stars.Where(star => star.SZ > 14)
                .ToList().Count;
            Assert.That(beforeOutsideSZRange, Is.Positive);
            var afterCount = stars.SelectOnSZ(minMag: null, maxMag: 14);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideSZRange));
            var afterOutsideSZRange = stars.Stars.Where(star => star.SZ > 14)
                .ToList().Count;
            Assert.That(afterOutsideSZRange, Is.Zero);
        }
    }
    
    [Test]
    public void Method_SelectOnSZ_Uncertainty_Test() {
        var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        var testUMag = 0.0105;  // 10.5 millimagnitudes.
        
        // Test cropping of SZ uncertainty to max limit:
        var beforeCount = stars.Count;
        var beforeUMagTooHigh = stars.Stars.Where(star => star.uSZ > testUMag)
            .ToList().Count;
        Assert.That(beforeUMagTooHigh, Is.Positive);
        var afterCount = stars.SelectOnSZ(maxUMag: testUMag);
        Assert.That(afterCount, Is.EqualTo(beforeCount - beforeUMagTooHigh));
        var afterUMagTooHigh = stars.Stars.Where(star => star.uSZ > testUMag)
            .ToList().Count;
        Assert.That(afterUMagTooHigh, Is.Zero);
    }
    
    [Test]
    public void Method_SelectOnRIColor_Test() {
        // Test cropping of R-I Color to range [0.2-0.3]:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideRIColorRange = stars.Stars.Where(star => star.RIColor is (< 0.2 or > 0.3))
                .ToList().Count;
            Assert.That(beforeOutsideRIColorRange, Is.Positive);
            var afterCount = stars.SelectOnRIColor(0.2, 0.3);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideRIColorRange));
            var afterOutsideRIColorRange = stars.Stars.Where(star => star.RIColor is (< 0.2 or > 0.3))
                .ToList().Count;
            Assert.That(afterOutsideRIColorRange, Is.Zero);
        }
        // Test cropping of R-I Color with min of 0.2 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideRIColorRange = stars.Stars.Where(star => star.RIColor < 0.2)
                .ToList().Count;
            Assert.That(beforeOutsideRIColorRange, Is.Positive);
            var afterCount = stars.SelectOnRIColor(minRIColor: 0.2, maxRIColor: null);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideRIColorRange));
            var afterOutsideRIColorRange = stars.Stars.Where(star => star.RIColor < 0.2)
                .ToList().Count;
            Assert.That(afterOutsideRIColorRange, Is.Zero);
        }
        // Test cropping of R-I Color with max of 0.3 only:
        {
            var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
            var beforeCount = stars.Count;
            var beforeOutsideRIColorRange = stars.Stars.Where(star => star.RIColor > 0.3)
                .ToList().Count;
            Assert.That(beforeOutsideRIColorRange, Is.Positive);
            var afterCount = stars.SelectOnRIColor(minRIColor: null, maxRIColor: 0.3);
            Assert.That(afterCount, Is.EqualTo(beforeCount - beforeOutsideRIColorRange));
            var afterOutsideRIColorRange = stars.Stars.Where(star => star.RIColor > 0.3)
                .ToList().Count;
            Assert.That(afterOutsideRIColorRange, Is.Zero);
        }
    }
    
    [Test]
    public void Method_ToString_Test() {
        var stars = new Atlas2Stars(cat, 23.75, 24.1, 33.7, 34.14, testDate);
        Assert.That(stars.ToString(), Is.EqualTo("Atlas2Stars: 639 star records."));
    }
} // class Atlas2StarsTests.
