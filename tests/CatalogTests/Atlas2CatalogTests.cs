using System.Diagnostics.CodeAnalysis;
using AstroLib.Catalog.Atlas2;
using NUnit.Framework;
using NUnit.Framework.Interfaces;

namespace AstroLibTests.CatalogTests;

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Atlas2CatalogTests {
    
    private string topPath;
    private Atlas2Catalog cat;

    [SetUp]
    public void Setup() {
        topPath = @"D:\Astro\Catalogs\ATLAS-refcat2";
        cat = new Atlas2Catalog(topPath);
    }

    [Test]
    public void Constructor_Tests() {
        Assert.IsInstanceOf(typeof(Atlas2Catalog), cat);
        Assert.That(cat.TopPath, Is.EqualTo(topPath));
        Assert.That(Atlas2Catalog.CatalogEpoch,
            Is.EqualTo(new DateTime(2015, 7, 1, 0, 0, 0)));
        Assert.That(Atlas2Catalog.RequiredFieldsPerDataLine, Is.EqualTo(44));
        Assert.That(Atlas2Catalog.SubdirectoryCodes.Count, Is.EqualTo(5));
        Assert.That(Atlas2Catalog.SubdirectoryCodes["mag-19-20"], Is.EqualTo("E"));
    }

    [Test]
    public void Property_Tests() {
        Assert.That(cat.TopPath, Is.EqualTo(topPath));
        Assert.That(cat.CatalogSubdirectoryPathsPresent.Count, Is.EqualTo(5));
        Assert.That(cat.CatalogSubdirectoryPathsPresent[4].StartsWith(topPath));
        Assert.That(cat.CatalogSubdirectoryPathsPresent[4].EndsWith("mag-19-20"));
        Assert.That(cat.CatalogSubdirectoryNamesPresent.Count, Is.EqualTo(5));
        Assert.That(cat.CatalogSubdirectoryNamesPresent[4], Is.EqualTo("mag-19-20"));
    }

    [Test]
    public void Method_Tests() {
        // .ValidateCatalogInstallation():
        var warnings = cat.ValidateCatalogInstallation();
        Assert.IsInstanceOf(typeof(List<string>), warnings);
        foreach (var w in warnings) {Console.WriteLine(w);}
        Assert.That(warnings.Count, Is.Zero);
        
        // .readLinesFromAtlasFile():
        var lines = cat.readLinesFromOneAtlasFile("mag-0-16", 149, 71);
        Assert.IsInstanceOf(typeof(List<string>), lines);
        Assert.That(lines.Count, Is.EqualTo(191));
        Assert.That(lines[98].Split(',').Count, 
            Is.EqualTo(Atlas2Catalog.RequiredFieldsPerDataLine));
        Assert.That(lines[98].StartsWith("14952694491,7137278481,58,2"));
        
        // .ToString():
        Assert.That(cat.ToString(), Is.EqualTo($"Atlas2Catalog object from {topPath}."));
    }
}