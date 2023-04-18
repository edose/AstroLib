using System.Diagnostics.CodeAnalysis;
using AstroLib.Core.Geometry;
using NUnit.Framework;

// ReSharper disable CompareOfFloatsByEqualityOperator
// ReSharper disable CommentTypo
// ReSharper disable InconsistentNaming
// ReSharper disable IdentifierTypo
// ReSharper disable InlineTemporaryVariable
// ReSharper disable JoinDeclarationAndInitializer

namespace AstroLibTests.GeometryTests; 

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class SphericalTests {

    [Test]
    public void MakeGoldenSpiral_Tests() {
        const int pointCount = 30;
        var pointList = Spherical.MakeGoldenSpiral(pointCount);
        Assert.That(pointList, Has.Count.EqualTo(pointCount));
        Assert.That(pointList, Has.Exactly(pointCount).Matches<SphCoord>(point => 
            point.Azimuth >= 0.0));
        Assert.That(pointList, Has.Exactly(pointCount).Matches<SphCoord>(point => 
            point.Azimuth < 2.0 * Math.PI));
        Assert.That(pointList, Has.Exactly(pointCount).Matches<SphCoord>(point => 
            point.PolarAngle >= 0.0));
        Assert.That(pointList, Has.Exactly(pointCount).Matches<SphCoord>(point => 
            point.PolarAngle <= Math.PI));
        Assert.That(pointList, Has.Exactly(pointCount).Matches<SphCoord>(point => 
            point.Radius == null));
        Assert.That(pointList[0].Azimuth, Is.EqualTo(5.083_204).Within(1E-6));
        Assert.That(pointList[0].PolarAngle, Is.EqualTo(0.258_922).Within(1E-6));
        Assert.That(pointList[27].Azimuth, Is.EqualTo(3.116_050).Within(1E-6));
        Assert.That(pointList[27].PolarAngle, Is.EqualTo(2.555_907).Within(1E-6));
    }
}

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class SphCoordTests {
    
    //TODO: Test that input values are within proper ranges. (HOW?)
    
    private readonly Point3D p1 = new Point3D(3, 4, 8);
    private readonly Vector3D v1 = new Vector3D(6.2, 8.9, 5.9);

    [Test]
    public void Constructor_Tests() {
        // Default constructor:
        // var sc1 = new SphCoord() {Azimuth = 1.3, PolarAngle = 1.4, Radius = 5.0};
        // Assert.That(sc1.Azimuth, Is.EqualTo(1.3));
        // Assert.That(sc1.PolarAngle, Is.EqualTo(1.4));
        // Assert.That(sc1.Radius, Is.EqualTo(5.0));
        // var sc2 = new SphCoord() {Azimuth = 2.3, PolarAngle = 2.4, Radius = null};
        // Assert.That(sc2.Azimuth, Is.EqualTo(2.3));
        // Assert.That(sc2.PolarAngle, Is.EqualTo(2.4));
        // Assert.That(sc2.Radius, Is.Null);

        // Constructor from 2 or 3 (with radius) doubles:
        var sc3 = new SphCoord(3.3, 3.4, 5.5);
        Assert.That(sc3.Azimuth, Is.EqualTo(3.3));
        Assert.That(sc3.PolarAngle, Is.EqualTo(3.4));
        Assert.That(sc3.Radius, Is.EqualTo(5.5));
        var sc4 = new SphCoord(3.3, 3.4, null);
        Assert.That(sc4.Azimuth, Is.EqualTo(3.3));
        Assert.That(sc4.PolarAngle, Is.EqualTo(3.4));
        Assert.That(sc4.Radius, Is.Null);
        var sc5 = new SphCoord(4.3, 4.4);
        Assert.That(sc5.Azimuth, Is.EqualTo(4.3));
        Assert.That(sc5.PolarAngle, Is.EqualTo(4.4));
        Assert.That(sc5.Radius, Is.Null);

        // Constructor from Point3D object:
        var sc_point = new SphCoord(p1);
        // TODO: write these tests.

        // Constructor from Vector3D object:
        // TODO: write these tests.

    }

    [Test]
    public void Property_Tests() {
        var sc1 = new SphCoord(3.3, 3.4, 5.5);
        Assert.That(sc1.Azimuth, Is.EqualTo(3.3));
        Assert.That(sc1.PolarAngle, Is.EqualTo(3.4));
        Assert.That(sc1.Radius, Is.EqualTo(5.5));
        
    }
    
    
    
}